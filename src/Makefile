.PHONY: all clean 
CXX=									g++
CXXFLAGS=							-g -Wall -Wno-unused-function -std=c++17 -O2
PROG=									../ms-index
INCLUDES-PATH?=~
INCLUDES=							-I$(INCLUDES-PATH)/include -L$(INCLUDES-PATH)/lib
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64

all:$(PROG)

clean:
	rm -f $(PROG)	

$(PROG): $(wildcard *.cpp *.c *.h)
	$(CXX) $(INCLUDES) $(CXXFLAGS) QSufSort.c main.cpp -o $@ $(LIBS) 

