# Makefile for compiling CN_search, precomputation, and Preproduct

# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -lgmp -pg

# Target executables
TARGETS = CN_search precomputation tab_job test small_tab

# Source files
SRCS = rollsieve.cpp CN_search.cpp precomputation.cpp Preproduct.cpp test_functions.cpp tab_job.cpp small_tabulation.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target (build all)
all: $(TARGETS)

# Rule for compiling CN_search
CN_search: CN_search.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling precomputation
precomputation: precomputation.o rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling tabulation job, shortened to tab_job
tab_job: tab_job.o Preproduct.o rollsieve.o
	$(CXX) $^  -o $@ $(CXXFLAGS)

# Rule for compiling test
test: test_functions.o Preproduct.o rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling small serial tabulation
small_tab:	small_tabulation.o Preproduct.o rollsieve.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Generic rule for compiling .cpp to .o
%.o:%.cpp
	$(CXX) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(OBJS) $(TARGETS)

# .PHONY to indicate these are not real files
.PHONY: all clean
