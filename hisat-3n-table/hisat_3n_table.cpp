/*
 * Single-threaded, fast I/O refactor of hisat-3n-table main driver.
 * Keeps original CLI (except threads are ignored), switches to fgets + streamless pipeline.
 */

 #include "position_3n_table.h"
 #include <getopt.h>
 #include <cstdio>
 #include <cstdlib>
 #include <cctype>
 #include <string>
 #include <iostream>
 #include <fstream>
 #include <chrono>
 #include <thread>
 
 using namespace std;
 
 // -------- Global variables --------
 long long int loadingBlockSize = 1000000LL;  // Default 1MB loading block size
 
 // -------- CLI state --------
 static string alignmentFileName;
 static bool   useStdin = true;
 static string refFileName;
 static string outputFileName;   // Let Positions decide stdout vs file
bool   uniqueOnly   = false;
bool   multipleOnly = false;
bool   CG_only      = false;
 // threads option is ignored in this refactor
char   convertFrom = '0';
char   convertTo   = '0';
char   convertFromComplement = 0;
char   convertToComplement   = 0;
bool   addedChrName   = false;
bool   removedChrName = false;
 
 // block size behavior is handled via meetNext; keep a default for clarity
 static const long long LOADING_BLOCK_SIZE = 1000000LL;
 
 // -------- utils --------
 static bool fileExist(const string& filename) {
     ifstream f(filename);
     return f.good();
 }
 
 enum {
     ARG_ADDED_CHRNAME = 256,
     ARG_REMOVED_CHRNAME
 };
 
 static const char *short_options = "a:r:o:b:umcp:h";
 static struct option long_options[] = {
     {"alignments",     required_argument, 0, 'a'},
     {"ref",            required_argument, 0, 'r'},
     {"output-name",    required_argument, 0, 'o'},
     {"base-change",    required_argument, 0, 'b'},
     {"unique-only",    no_argument,       0, 'u'},
     {"multiple-only",  no_argument,       0, 'm'},
     {"CG-only",        no_argument,       0, 'c'},
     {"threads",        required_argument, 0, 'p'}, // ignored in single-thread mode
     {"added-chrname",  no_argument,       0, ARG_ADDED_CHRNAME},
     {"removed-chrname",no_argument,       0, ARG_REMOVED_CHRNAME},
     {"help",           no_argument,       0, 'h'},
     {0, 0, 0, 0}
 };
 
 static void printHelp(ostream& out) {
     out << "hisat-3n-table (single-thread, fast I/O)\n"
         << "Usage:\n"
         << "  hisat-3n-table [options]* -a <alignment.sam|-> -r <ref.fa> -b <X,Y> [-o <out.tsv>]\n"
         << "  <alignment.sam|->  SORTED SAM filename, or '-' for stdin.\n"
         << "  <ref.fa>           FASTA reference.\n"
         << "  <X,Y>              --base-change (e.g. C,T)\n"
         << "Options:\n"
         << "  -u/--unique-only       count only unique-mapped bases\n"
         << "  -m/--multiple-only     count only multi-mapped bases (mutually exclusive with -u)\n"
         << "  -c/--CG-only           limit to CG sites (forces --base-change=C,T)\n"
         << "  --added-chrname        set if alignment used --add-chrname\n"
         << "  --removed-chrname      set if alignment used --remove-chrname\n"
         << "  -o/--output-name FILE  output TSV file (default: stdout)\n"
         << "  -p/--threads <int>     (ignored in this build)\n"
         << "  -h/--help              show this help\n";
 }
 
 static void parseOption(int next_option, const char *optarg) {
     switch (next_option) {
         case 'a': {
             alignmentFileName = optarg;
             useStdin = (alignmentFileName == "-");
             if (!useStdin && !fileExist(alignmentFileName)) {
                 cerr << "The alignment file does not exist.\n";
                 throw 1;
             }
             break;
         }
         case 'r': {
             refFileName = optarg;
             if (!fileExist(refFileName)) {
                 cerr << "Reference (FASTA) file does not exist.\n";
                 throw 1;
             }
             break;
         }
         case 'o': {
             outputFileName = optarg; // Let Positions decide to open this path or stdout
             break;
         }
         case 'b': {
             string arg = optarg ? string(optarg) : string();
             if (arg.size() != 3 || arg[1] != ',') {
                 cerr << "Error: expected 'X,Y' for --base-change (e.g. C,T), got " << arg << "\n";
                 throw 1;
             }
             convertFrom = toupper(arg.front());
             convertTo   = toupper(arg.back());
             break;
         }
         case 'u': uniqueOnly = true; break;
         case 'm': multipleOnly = true; break;
         case 'c': CG_only = true; break;
         case 'p': /* ignore threads */ break;
         case ARG_ADDED_CHRNAME:   addedChrName = true; break;
         case ARG_REMOVED_CHRNAME: removedChrName = true; break;
         case 'h': printHelp(cerr); throw 0;
         default:  printHelp(cerr); throw 1;
     }
 }
 
 static void parseOptions(int argc, const char **argv) {
     int option_index = 0;
     int next_option;
     while (true) {
         next_option = getopt_long(argc, const_cast<char **>(argv), short_options,
                                   long_options, &option_index);
         if (next_option == -1) break;
         parseOption(next_option, optarg);
     }
 
     if (refFileName.empty() || alignmentFileName.empty()) {
         cerr << "No reference or SAM file specified!\n";
         printHelp(cerr);
         throw 1;
     }
 
     if (CG_only) {
         if (convertFrom != 'C' || convertTo != 'T') {
             cerr << "Warning: CG-only mode enforces --base-change=C,T\n";
             convertFrom = 'C';
             convertTo   = 'T';
         }
     }
 
     if (convertFrom == '0' || convertTo == '0') {
         cerr << "Error: --base-change is required.\n";
         throw 1;
     }
 
     if (removedChrName && addedChrName) {
         cerr << "Error: --removed-chrname and --added-chrname cannot be used together.\n";
         throw 1;
     }
     if (uniqueOnly && multipleOnly) {
         cerr << "Error: --unique-only and --multiple-only are mutually exclusive.\n";
         throw 1;
     }
 
     // complements (assumes asc2dnacomp[] available from position_3n_table.h)
     convertFromComplement = asc2dnacomp[convertFrom];
     convertToComplement   = asc2dnacomp[convertTo];
 }
 
 // Extract RNAME (field 3) and POS (field 4) from a SAM line; return false if unmapped.
 static bool getSAMChromosomePos(const string& line, string& chr, long long& pos) {
     int start = 0, end = 0, count = 0;
     while ((end = line.find('\t', start)) != (int)string::npos) {
         if (count == 2) {
             chr = line.substr(start, end - start);
         } else if (count == 3) {
             pos = stoll(line.substr(start, end - start));
             if (chr == "*") return false;
             return true;
         }
         start = end + 1;
         ++count;
     }
     return false;
 }
 
 static int hisat_3n_table() {
     // ---- Fast I/O setup ----
     ios::sync_with_stdio(false);
 
     // Open input stream for fgets
     FILE* fin = stdin;
     if (!useStdin) {
         fin = fopen(alignmentFileName.c_str(), "r");
         if (!fin) {
             perror("fopen");
             return 1;
         }
     }
 
         // ---- Construct Positions (single-thread mode) ----
    Positions positions(refFileName);
    
    // Start output thread
    thread outputThread(&Positions::outputFunction, &positions, outputFileName);
    outputThread.detach();
 
     // ---- Stream pipeline ----
     static char buff[1000007]; // large line buffer for fgets
     string samChromosome;
     long long samPos = 0;
     long long reloadPos = 0;       // when samPos > reloadPos => loadMore
     long long lastPos   = 0;       // verify SAM is coordinate-sorted within a chromosome
     const int INF = 0x3f3f3f3f;
 
     bool firstChromLoaded = false;
 
     while (true) {
         if (!fgets(buff, sizeof(buff), fin)) break;
         string line(buff);
 
         if (line.empty() || line[0] == '@') continue;
 
         if (!getSAMChromosomePos(line, samChromosome, samPos)) continue;
 
         // First time / chromosome changed
         if (!firstChromLoaded || samChromosome != positions.chromosome) {
             // flush everything before switching
             positions.startOutput(true);
             int meetNext = 0;
             positions.loadNewChromosome(samChromosome, meetNext);
             reloadPos = meetNext ? INF : LOADING_BLOCK_SIZE;
             lastPos   = 0;
             firstChromLoaded = true;
         }
 
         // Need to extend reference window
         while (samPos > reloadPos) {
             positions.startOutput(); // flush block output if any
             int meetNext = 0;
             positions.loadMore(meetNext);
             reloadPos += meetNext ? INF : LOADING_BLOCK_SIZE;
         }
 
         if (lastPos > samPos) {
             cerr << "The input alignment file is not sorted; please provide a coordinate-sorted SAM.\n";
             if (fin && fin != stdin) fclose(fin);
             return 1;
         }
 
         // Synchronous append (no queues, no workers)
         positions.appendSync(line);
         lastPos = samPos;
     }
 
         // Final flush and stop working
    positions.startOutput(true);
    positions.working = false;
    
    // Wait a bit for output to complete
    this_thread::sleep_for(chrono::milliseconds(100));

    if (fin && fin != stdin) fclose(fin);
    return 0;
 }
 
 int main(int argc, const char** argv) {
     try {
         parseOptions(argc, argv);
         return hisat_3n_table();
     } catch (const std::exception& e) {
         cerr << "Error: exception: '" << e.what() << "'\nCommand: ";
         for (int i = 0; i < argc; ++i) cerr << argv[i] << ' ';
         cerr << '\n';
         return 1;
     } catch (int e) {
         if (e != 0) {
             cerr << "Error: internal HISAT-3N exception (#" << e << ")\nCommand: ";
             for (int i = 0; i < argc; ++i) cerr << argv[i] << ' ';
             cerr << '\n';
         }
         return e;
     }
 }
 