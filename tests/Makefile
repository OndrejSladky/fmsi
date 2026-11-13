.PHONY: all clean test integration-test

CXX = g++
CXXFLAGS= -g -Wall -Wno-unused-function -std=c++17 -O2
INCLUDES-PATH?=../src
INCLUDES=							-I$(INCLUDES-PATH)/include -L$(INCLUDES-PATH)/lib
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64
GTEST= googletest/googletest
SRC= ../src
LIBSAIS= $(SRC)/libsais/libsais.o $(SRC)/libsais/libsais64.o

CC?= cc
include $(SRC)/openmp.mk
ifneq ($(OPENMP_CFLAGS),)
CXXFLAGS+= -DLIBSAIS_OPENMP
endif


all: test integration-test verification-test

test: fmsi-test
	./fmsi-test

integration-test: ./integration_test.sh
	$(MAKE) -C ..
	./integration_test.sh

verification-test: ./verify.py
	$(MAKE) -C ..
	python3 verify.py data/GCF_009858895.2_ASM985889v3_genomic.fna.ms.k13 --num_queries=100 --k=13 --streaming_length=50
	python3 verify.py data/GCF_test2_ones.fna.ms.k31 --num_queries=100 --k=31 --streaming_length=60 --optimize=True
	python3 verify.py data/GCF_009858895.2_ASM985889v3_genomic.fna.ms.k13 --k=13 --hashes=True

fmsi-test: $(wildcard *.cpp *.h *.hpp) gtest-all.o $(wildcard $(SRC)/*.cpp $(SRC)/*.h) $(SRC)/QSufSort.c $(LIBSAIS) $(OPENMP_SETTING)
	$(CXX) $(CXXFLAGS) -I $(GTEST)/include $(INCLUDES) main.cpp $(SRC)/QSufSort.c $(LIBSAIS) gtest-all.o -pthread -o $@ $(LIBS) $(OPENMP_LIBS)

$(LIBSAIS):
	$(MAKE) -C $(SRC) $(LIBSAIS:$(SRC)/%=%)

gtest-all.o: $(GTEST)/src/gtest-all.cc
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f fmsi-test
	rm -f gtest-all.o
