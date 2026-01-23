# Makefile for compiling CN_search, precomputation, and Preproduct

# Compiler and flags.  Note that -I. tells compiler to start in current directory for relative paths
CXX = g++
CXXFLAGS = -O3 -lgmp -lgmpxx -I.
MPI = mpic++
PATHS = -I/usr/local/include -L/usr/local/lib

# Target executables
TARGETS = CN_query precomputation tab_job test small_tab parallel_large parallel_small

# Source files
SRCS = IncrementalSieve/rollsieve.cpp CN_query.cpp TabulationFiles/precomputation.cpp Preproduct.cpp TabulationFiles/test_functions.cpp TabulationFiles/tab_job.cpp TabulationFiles/small_tabulation.cpp TabulationFiles/parallel_small_P.cpp TabulationFiles/parallel_large_P.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target (build all)
all: $(TARGETS)

# Generic rule for compiling .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for compiling CN_query
CN_query: CN_query.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling precomputation
precomputation: TabulationFiles/precomputation.o IncrementalSieve/rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling tabulation job, shortened to tab_job
tab_job: TabulationFiles/tab_job.o Preproduct.o IncrementalSieve/rollsieve.o
	$(CXX) $^  -o $@ $(CXXFLAGS)

# Rule for compiling test
test: TabulationFiles/test_functions.o Preproduct.o IncrementalSieve/rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling small serial tabulation
small_tab: TabulationFiles/small_tabulation.o Preproduct.o IncrementalSieve/rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling parallel_large_P
parallel_large: IncrementalSieve/rollsieve.o Preproduct.o
	mpic++ -O3 -lgmp -lgmpxx -I. -c TabulationFiles/parallel_large_P.cpp -o TabulationFiles/parallel_large_P.o

# Rule for compiling parallel_small_P
parallel_small:  IncrementalSieve/rollsieve.o Preproduct.o
	$(MPI) $(PATHS) $^ -o $@ $(CXXFLAGS)

# Clean up object files and executables
clean:
	rm -f $(OBJS) $(TARGETS)

# .PHONY to indicate these are not real files
.PHONY: all clean
