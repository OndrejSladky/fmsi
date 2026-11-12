.PHONY: all clean 
CXX=									g++
ARCH?=								native
CXXFLAGS=							-g -Wall -Wno-unused-function -std=c++17 -O3 -march=$(ARCH) -pthread
CC?=									cc
CFLAGS=								-g -O2
PROG=									../fmsi
INCLUDES-PATH?=.
INCLUDES=							-I$(INCLUDES-PATH)/include -L$(INCLUDES-PATH)/lib 
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64
LIBSAIS=							libsais/libsais.o libsais/libsais64.o

all:$(PROG)


include/sdsl/suffix_arrays.hpp: sdsl-lite/build/Makefile
	$(MAKE) -C sdsl-lite/build
	$(MAKE) -C sdsl-lite/build install

sdsl-lite/build/Makefile:
	cd sdsl-lite/build && cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=../.. ..


$(PROG): $(wildcard *.cpp *.c *.h) $(LIBSAIS) ./include/sdsl/suffix_arrays.hpp version.h
	./create-version.sh
	$(CXX) $(INCLUDES) $(CXXFLAGS) QSufSort.c main.cpp $(LIBSAIS) -o $@ $(LIBS)

libsais/%.o: libsais/%.c libsais/libsais.h libsais/libsais64.h
	$(CC) $(CFLAGS) -c $< -o $@

version.h: version
	./create-version.sh


clean:
	rm -f $(PROG)
	rm -f version.h
	rm -f $(LIBSAIS)
	$(MAKE) -C sdsl-lite/build/ clean
