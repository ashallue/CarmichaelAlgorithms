# Makefile for compiling CN_search, precomputation, and Preproduct

# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -lgmp -lgmpxx

# Target executables
TARGETS = CN_search precomputation Preproduct

# Source files
SRCS = CN_search.cpp precomputation.cpp Preproduct.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target (build all)
all: $(TARGETS)

# Rule for compiling CN_search
CN_search: CN_search.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling precomputation
precomputation: precomputation.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Rule for compiling Preproduct
Preproduct: Preproduct.o
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Generic rule for compiling .cpp to .o
%.o: %.cpp
	$(CXX) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(OBJS) $(TARGETS)

# .PHONY to indicate these are not real files
.PHONY: all clean
