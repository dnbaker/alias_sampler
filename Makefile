CXX?=g++

STD?=c++14
OPT?=-O3

CXXFLAGS = -std=$(STD) $(OPT) -fno-inline

all: test
%: %.cpp $(wildcard *.h)
	$(CXX) $(CXXFLAGS) $< -o $@
