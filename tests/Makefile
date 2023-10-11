.PHONY: all clean test

CXX = g++
CXXFLAGS= -g -Wall -Wno-unused-function -std=c++17 -O2
GTEST= googletest/googletest
SRC= ../src

all: test

test: ms-index-test
	./ms-index-test

ms-index-test: $(wildcard *.cpp *.h *.hpp) gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include main.cpp gtest-all.o -pthread -o $@

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f ms-index-test
	rm -f gtest-all.o
