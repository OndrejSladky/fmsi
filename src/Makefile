.PHONY: all clean all-ci-linux all-ci-macos
CXX=									g++
CXXFLAGS=							-g -Wall -Wno-unused-function -std=c++17 -O2
PROG=									../ms-index
INCLUDES=							-I./include -L./lib
INCLUDES-CI-LINUX=		-I/home/runner/include -L/home/runner/lib
INCLUDES-CI-MACOS=		-I/Users/runner/include -L/Users/runner/lib
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64

all:$(PROG) ./include/sdsl/suffix_arrays.hpp

all-ci-linux:
	$(CXX) $(INCLUDES-CI-LINUX) $(CXXFLAGS)  QSufSort.c main.cpp -o $@ $(LIBS)

all-ci-macos:
	$(CXX) $(INCLUDES-CI-MACOS) $(CXXFLAGS)  QSufSort.c main.cpp -o $@ $(LIBS)

include/sdsl/suffix_arrays.hpp: sdsl-lite/build/Makefile
	$(MAKE) -C sdsl-lite/build
	$(MAKE) -C sdsl-lite/build install

sdsl-lite/build/Makefile:
	cd sdsl-lite/build && cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=../.. ..

clean:
	rm -f $(PROG)
	$(MAKE) -C sdsl-lite/build/ clean


$(PROG): $(wildcard *.cpp *.c *.h) ./include/sdsl/suffix_arrays.hpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) QSufSort.c main.cpp -o $@ $(LIBS)

