.PHONY: all clean 
CXX=									g++
CXXFLAGS=							-g -Wall -Wno-unused-function -std=c++17 -O2
PROG=									../fmsi
INCLUDES-PATH?=.
INCLUDES=							-I$(INCLUDES-PATH)/include -L$(INCLUDES-PATH)/lib 
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64

all:$(PROG)


include/sdsl/suffix_arrays.hpp: sdsl-lite/build/Makefile
	$(MAKE) -C sdsl-lite/build
	$(MAKE) -C sdsl-lite/build install

sdsl-lite/build/Makefile:
	cd sdsl-lite/build && cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=../.. ..


$(PROG): $(wildcard *.cpp *.c *.h) ./include/sdsl/suffix_arrays.hpp version.h
	./create-version.sh
	$(CXX) $(INCLUDES) $(CXXFLAGS) QSufSort.c main.cpp -o $@ $(LIBS)

version.h: version
	./create-version.sh


clean:
	rm -f $(PROG)
	rm -f version.h
	$(MAKE) -C sdsl-lite/build/ clean
