#include <htslib/sam.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct Rec {
    bam1_t* b = nullptr;
    int32_t tid = -1;
    int64_t pos = -1;
    Rec() = default;
    explicit Rec(const bam1_t* src) {
        b = bam_dup1(src);
        tid = b->core.tid;
        pos = b->core.pos;
    }
    // 按位点排序
    bool operator<(const Rec& other) const {
        return pos < other.pos;
    }
};

static void die(const string& msg) {
    cerr << "[ERROR] " << msg << "\n";
    exit(1);
}

static bool headers_compatible(const bam_hdr_t* a, const bam_hdr_t* b) {
    if (a->n_targets != b->n_targets) return false;
    for (int i = 0; i < a->n_targets; ++i) {
        if (strcmp(a->target_name[i], b->target_name[i]) != 0) return false;
        if (a->target_len[i] != b->target_len[i]) return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0]
             << " <input_prefix> <input_num> <output_num> [--out BAM|SAM]\n";
        return 1;
    }
    const string prefix = argv[1];
    const int input_num = stoi(argv[2]);
    const int output_num = stoi(argv[3]);
    string out_fmt = "BAM"; // or SAM

    for (int i = 4; i < argc; ++i) {
        string s = argv[i];
        if (s == "--out" && i+1 < argc) {
            out_fmt = argv[++i];
            if (out_fmt != "BAM" && out_fmt != "SAM") {
                die("Invalid --out (BAM|SAM)");
            }
        } else {
            die("Unknown argument: " + s);
        }
    }
    if (input_num <= 0 || output_num <= 0) die("input_num and output_num must be > 0");

    // ---------- 读取第一个文件，建立“主”header ----------
    char path[1024];
    snprintf(path, sizeof(path), "%s.%d.bam", prefix.c_str(), 0);
    htsFile* fp0 = sam_open(path, "rb");
    if (!fp0) die("Failed to open " + string(path));
    bam_hdr_t* master_hdr = sam_hdr_read(fp0);
    if (!master_hdr) die("Failed to read header from " + string(path));
    sam_close(fp0);

    const int n_targets = master_hdr->n_targets;
    cerr << "[INFO] n_targets = " << n_targets << "\n";

    // ---------- 每条染色体一个桶（暂存全部 reads） ----------
    vector<vector<Rec>> chr_vec(n_targets);

    // ---------- 读入所有输入文件 ----------
    cerr << "[INFO] Loading inputs into memory...\n";
    for (int i = 0; i < input_num; ++i) {
        snprintf(path, sizeof(path), "%s.%d.bam", prefix.c_str(), i);
        htsFile* fp = sam_open(path, "rb");
        if (!fp) die("Failed to open " + string(path));
        bam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) die("Failed to read header from " + string(path));

        if (!headers_compatible(master_hdr, hdr)) {
            die("Headers are not compatible across inputs (different target names/lengths).");
        }

        bam1_t* b = bam_init1();
        while (sam_read1(fp, hdr, b) >= 0) {
            if (b->core.tid < 0 || b->core.tid >= n_targets) continue; // unmapped等
            chr_vec[b->core.tid].emplace_back(b);
        }
        bam_destroy1(b);
        sam_hdr_destroy(hdr);
        sam_close(fp);
    }

    // ---------- 统计 reads/染色体 & 过滤空染色体 ----------
    struct ChrStat { int tid; uint64_t reads; };
    vector<ChrStat> chr_stats;
    chr_stats.reserve(n_targets);
    for (int tid = 0; tid < n_targets; ++tid) {
        if (!chr_vec[tid].empty()) {
            chr_stats.push_back({tid, static_cast<uint64_t>(chr_vec[tid].size())});
        }
    }
    // 按 reads 降序
    sort(chr_stats.begin(), chr_stats.end(), [](const ChrStat& a, const ChrStat& b){
        return a.reads > b.reads;
    });

    cerr << "[INFO] Chromosome distribution (non-empty only):\n";
    for (auto &cs : chr_stats) {
        cerr << "  " << master_hdr->target_name[cs.tid] << ": " << cs.reads << " reads\n";
    }

    // ---------- round-robin 分配到 output_num 个“桶” ----------
    vector<vector<int>> bucket_chrs(output_num); // 每个输出桶负责哪些 tid
    for (size_t i = 0; i < chr_stats.size(); ++i) {
        int b = static_cast<int>(i % output_num);
        bucket_chrs[b].push_back(chr_stats[i].tid);
    }

    cerr << "[INFO] Mapping (bucket -> chromosome list):\n";
    for (int b = 0; b < output_num; ++b) {
        uint64_t sum = 0;
        for (int tid : bucket_chrs[b]) sum += chr_vec[tid].size();
        cerr << "  output_" << b << " (" << sum << " reads):";
        for (int tid : bucket_chrs[b]) cerr << " " << master_hdr->target_name[tid];
        cerr << "\n";
    }

    // ---------- 打开输出文件 ----------
    struct OutHandle {
        htsFile* fp = nullptr;
        bam_hdr_t* hdr = nullptr; // 独立一份 header（避免多线程共享写时隐患）
    };
    vector<OutHandle> outs(output_num);

    for (int b = 0; b < output_num; ++b) {
        string ext = (out_fmt == "BAM") ? "bam" : "sam";
        string out_path = prefix + ".sorted.split." + to_string(b) + "." + ext;
        const char* mode = (out_fmt == "BAM") ? "wb" : "w"; // BAM 二进制 / SAM 文本

        outs[b].fp = sam_open(out_path.c_str(), mode);
        if (!outs[b].fp) die("Failed to open output " + out_path);
        outs[b].hdr = bam_hdr_dup(master_hdr);
        if (!outs[b].hdr) die("bam_hdr_dup failed");
        if (sam_hdr_write(outs[b].fp, outs[b].hdr) < 0) {
            die("sam_hdr_write failed for " + out_path);
        }
    }

    // ---------- 并行写：每个线程负责一个输出桶，桶内串行按染色体写 ----------
    cerr << "[INFO] Writing outputs in parallel...\n";
    #pragma omp parallel for schedule(dynamic,1)
    for (int b = 0; b < output_num; ++b) {
        // 本线程只写 outs[b]
        auto* fp = outs[b].fp;
        auto* hdr = outs[b].hdr;

        for (int tid : bucket_chrs[b]) {
            auto& vec = chr_vec[tid];
            // 对该染色体的 reads 按位置排序
            sort(vec.begin(), vec.end()); // Rec::operator<

            for (auto& r : vec) {
                if (sam_write1(fp, hdr, r.b) < 0) {
                    #pragma omp critical
                    {
                        cerr << "[ERROR] sam_write1 failed on bucket " << b
                             << " chr " << hdr->target_name[tid] << "\n";
                    }
                    // 不中断其他线程，可继续尝试
                }
            }
            #pragma omp critical
            {
                cerr << "[INFO] bucket " << b << " wrote chr "
                     << hdr->target_name[tid] << " (" << vec.size() << ")\n";
            }
        }
    }

    // ---------- 关闭输出，释放资源 ----------
    for (int b = 0; b < output_num; ++b) {
        if (outs[b].fp) sam_close(outs[b].fp);
        if (outs[b].hdr) sam_hdr_destroy(outs[b].hdr);
    }

    // 释放内存中所有 bam1_t
    for (int tid = 0; tid < n_targets; ++tid) {
        for (auto& r : chr_vec[tid]) {
            if (r.b) bam_destroy1(r.b);
        }
    }
    sam_hdr_destroy(master_hdr);

    cerr << "[INFO] Done.\n";
    return 0;
}
