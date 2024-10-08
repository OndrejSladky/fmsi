.PHONY: all clean test integration-test

CXX = g++
CXXFLAGS= -g -Wall -Wno-unused-function -std=c++17 -O2
INCLUDES-PATH?=../src
INCLUDES=							-I$(INCLUDES-PATH)/include -L$(INCLUDES-PATH)/lib
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64
GTEST= googletest/googletest
SRC= ../src


all: test integration-test verification-test

test: fmsi-test
	./fmsi-test

integration-test: ./integration_test.sh
	$(MAKE) -C ..
	./integration_test.sh

verification-test: ./verify.py
	$(MAKE) -C ..
	python3 verify.py data/GCF_009858895.2_ASM985889v3_genomic.fna.ms.k13 --num_queries=1000 --k=13

fmsi-test: $(wildcard *.cpp *.h *.hpp) gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(SRC)/QSufSort.c
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include  $(INCLUDES) main.cpp $(SRC)/QSufSort.c gtest-all.o -pthread -o $@ $(LIBS)

gtest-all.o: $(GTEST)/src/gtest-all.cc
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f fmsi-test
	rm -f gtest-all.o
