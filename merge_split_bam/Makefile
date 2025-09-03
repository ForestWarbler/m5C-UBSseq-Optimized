HTSLIB_DIR := /GLOBALFS/sysu_hpcscc_2/zhj/projects/ASC25_simulation/rna_2/samtools-1.21/htslib-1.21/build
CXX      := g++
CXXFLAGS := -O2 -std=c++17 -fopenmp -I$(HTSLIB_DIR)/include
LDFLAGS  := -L$(HTSLIB_DIR)/lib -lhts

TARGET   := merge_split_bam
SRC      := merge_split_bam_htslib_omp.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o
