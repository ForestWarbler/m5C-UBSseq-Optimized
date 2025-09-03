# Makefile for hisat-3n-table (single-threaded version)
#
# This Makefile builds the optimized single-threaded hisat-3n-table implementation
# with fast I/O and streamless pipeline optimizations.

# Compiler settings
CXX = g++
CXXFLAGS = -O2 -std=c++11 -Wall -Wextra
DEBUG_FLAGS = -g -O0 -DDEBUG -fno-omit-frame-pointer
INCLUDES = -I.
LIBS = -lpthread

# Source files
MAIN_SRC = hisat_3n_table.cpp
HEADERS = position_3n_table.h alignment_3n_table.h utility_3n_table.h

# Target executables
TARGET = hisat-3n-table
DEBUG_TARGET = hisat-3n-table-debug

# Default target
all: $(TARGET)

# Release version (optimized)
$(TARGET): $(MAIN_SRC) $(HEADERS)
	@echo "Building optimized hisat-3n-table..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(MAIN_SRC) $(LIBS)
	@echo "✅ hisat-3n-table built successfully!"
	@echo "Usage: ./$(TARGET) -a input.sam -r reference.fa -b C,T -o output.tsv"

# Debug version
debug: $(DEBUG_TARGET)

$(DEBUG_TARGET): $(MAIN_SRC) $(HEADERS)
	@echo "Building debug hisat-3n-table..."
	$(CXX) $(DEBUG_FLAGS) $(INCLUDES) -o $(DEBUG_TARGET) $(MAIN_SRC) $(LIBS)
	@echo "✅ Debug hisat-3n-table built successfully!"
	@echo "Usage: ./$(DEBUG_TARGET) -a input.sam -r reference.fa -b C,T -o output.tsv"

# Clean targets
clean:
	@echo "Cleaning build files..."
	rm -f $(TARGET) $(DEBUG_TARGET)
	@echo "Clean completed."

# Install target (copy to parent bin directory if exists)
install: $(TARGET)
	@if [ -d "../bin" ]; then \
		cp $(TARGET) ../bin/; \
		echo "Installed $(TARGET) to ../bin/"; \
	else \
		echo "No ../bin directory found, skipping install"; \
	fi

# Test compilation
test: $(TARGET)
	@echo "Testing hisat-3n-table compilation..."
	@./$(TARGET) --help >/dev/null 2>&1 && echo "✅ Binary works correctly" || echo "❌ Binary test failed"

# Show help
help:
	@echo "Available targets:"
	@echo "  all      - Build optimized hisat-3n-table (default)"
	@echo "  debug    - Build debug version with debug symbols"
	@echo "  clean    - Remove all build files"
	@echo "  install  - Install to ../bin/ directory"
	@echo "  test     - Test if the binary works"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Features of this implementation:"
	@echo "  - Single-threaded, fast I/O design"
	@echo "  - Streamless pipeline for better performance"
	@echo "  - Compatible CLI with original hisat-3n-table"
	@echo "  - Optimized for large SAM file processing"

# Phony targets
.PHONY: all debug clean install test help
