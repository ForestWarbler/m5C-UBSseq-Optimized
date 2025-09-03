# m5C-UBSseq-Optimized

An optimized bioinformatics pipeline for processing 5-methylcytosine (m5C) data from UBS-seq (Uracil-Bisulfite-Sequencing) experiments. This pipeline implements read-level parallelization and intelligent chromosome splitting for high-performance analysis of RNA methylation sequencing data.

## Features

- **Read-level Parallelization**: Splits input FASTQ files into chunks for parallel processing
- **Intelligent Chromosome Splitting**: Uses custom C++ tools for efficient BAM file splitting
- **Contamination Filtering**: Two-stage mapping (contamination → genome) for clean data
- **UMI-based Deduplication**: Support for UMI-based duplicate removal using UMICollapse
- **Quality Control**: Comprehensive quality filtering and reporting
- **Methylation Site Detection**: Specialized calling pipeline for m5C sites

## Pipeline Stages

### Stage 1: Read Processing and Mapping (Snakefile-stage1)
1. **Read Trimming**: Quality trimming and adapter removal using `cutseq`
2. **Read Splitting**: Parallel processing by splitting reads into chunks
3. **Contamination Mapping**: Initial mapping against contamination reference
4. **Genome Mapping**: Mapping unmapped reads against genome reference
5. **Deduplication**: UMI-based or mark-duplicate based deduplication
6. **Methylation Calling**: Site-specific methylation detection

### Stage 2: Site Analysis and Filtering (stage2.sh)
1. **Group Analysis**: Combine results from multiple samples
2. **Site Selection**: Identify potential methylation sites
3. **Background Filtering**: Filter sites based on background noise
4. **Final Output**: Generate filtered methylation site lists

## Dependencies

### Core Tools (included as submodules)
- **HISAT-3N**: Modified HISAT2 for bisulfite sequencing alignment
- **Samtools**: SAM/BAM file manipulation
- **UMICollapse**: UMI-based duplicate removal
- **SRA-tools**: SRA data access utilities
- **libdeflate**: Compression library

### Python Environment
The pipeline requires a conda environment with bioinformatics packages. See `environment.yml` for details.

### Custom Tools
- **merge_split_bam**: Custom C++ tool for intelligent BAM file splitting with OpenMP parallelization
- **hisat-3n-table**: Custom optimized single-threaded implementation with fast I/O and streamless pipeline

## Installation

1. **Clone the repository with submodules**:
   ```bash
   git clone --recursive https://github.com/your-repo/m5C-UBSseq-Optimized.git
   cd m5C-UBSseq-Optimized
   ```

2. **Set up the conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate myenv
   ```

3. **Build the submodules**:
   ```bash
   # Build HISAT-3N
   cd hisat-3n
   make -j$(nproc)
   cd ..
   
   # Build Samtools
   cd samtools
   ./configure --prefix=$(pwd)
   make -j$(nproc)
   cd ..
   
   # Build UMICollapse
   cd UMICollapse
   # Follow UMICollapse build instructions
   cd ..
   
   # Build custom merge_split_bam tool
   cd merge_split_bam
   make -j$(nproc)
   cd ..
   
   # Build custom hisat-3n-table tool
   cd hisat-3n-table
   make
   cd ..
   ```

4. **Build SRA-tools and libdeflate** (if needed):
   ```bash
   cd sra-tools
   ./configure --prefix=$(pwd)
   make -j$(nproc)
   cd ..
   
   cd libdeflate
   make -j$(nproc)
   cd ..
   ```

## Configuration

Create a `config.yaml` file with the following structure:

```yaml
# tools and scripts
path:
  samtools: path/to/samtools 
  hisat3n: path/to/hisat-3n
  hisat3ntable: path/to/hisat-3n-table
  umicollapse: path/to/umicollapse.jar
  join_pileup.py: path/to/join_pileup.py
  group_pileup.py: path/to/group_pileup.py
  select_sites.py: path/to/select_sites.py
  filter_sites.py: path/to/filter_sites.py
  merge_split_bam.py: path/to/merge_split_bam.py

# global config
#
# prepare genes index
# premap to rRNA, tRNA and other small RNA
# If study virus, then also premap to virus genome
# customized_genes:
#   - a.fa
#   - b.fa
library: INLINE
# makedup: false

# reference genome and index
reference:
  contamination:
    fa: path/to/contamination.fa
    hisat3n: path/to/contamination
  genes:
    fa: path/to/Homo_sapiens.GRCh38.sncRNA.fa
    hisat3n: path/to/Homo_sapiens.GRCh38.sncRNA
  genome:
    fa: path/to/Homo_sapiens.GRCh38.genome.fa
    hisat3n: path/to//Homo_sapiens.GRCh38.genome


# Sample name should be indentical and listed in the 2nd level of the yaml file
# Each sample will be analysis seperately, but
# samples sharing the same group id will be regared as biological replicates and combined in the comparing step
samples:
  SRRXXXXXXXX:
    data:
      - R1: path/to/SRRXXXXXXXX_1.fastq.gz
    group: TEST1
```

## Usage

### Stage 1: Primary Processing
```bash
# Run the main pipeline with read-level parallelization
snakemake -s Snakefile-stage1 --cores 64 --use-conda
```

### Stage 2: Site Analysis
```bash
# Configure stage2.sh with your paths
nano stage2.sh

# Run stage 2 analysis
bash stage2.sh
```

## Performance Optimizations

### Read-level Parallelization
- Input FASTQ files are split into 20 chunks by default
- Each chunk is processed independently through the entire mapping pipeline
- Results are merged back together for downstream analysis

### Intelligent Chromosome Splitting
- Custom C++ tool (`merge_split_bam`) implements smart BAM file splitting
- Uses OpenMP for parallel processing
- Optimizes memory usage and I/O operations

### Custom hisat-3n-table Implementation
- **Single-threaded Optimization**: Designed for reliable single-threaded performance
- **Fast I/O Pipeline**: Uses `fgets` instead of stream-based processing for better performance
- **Streamless Processing**: Eliminates memory overhead from streaming operations
- **Compatible Interface**: Maintains full compatibility with original hisat-3n-table CLI
- **Memory Efficient**: Optimized for processing large SAM files with minimal memory footprint

### Resource Management
- Configurable thread allocation per rule
- Temporary file cleanup to manage disk space
- Memory-efficient processing for large datasets
