.PHONY: all clean
CC=					g++
CFLAGS=			-g -Wall -Wno-unused-function -g3
PROG=				../ms-index
INCLUDES=		-I/usr/local/include -L/usr/local/lib
LIBS=				-lz -lsdsl -ldivsufsort -ldivsufsort64

all:$(PROG)

clean:
	rm -f $(PROG)	

$(PROG): $(wildcard *.cpp *.c *.h)
	$(CC) $(INCLUDES) $(CFLAGS) $(DFLAGS) $(AOBJS2) QSufSort.c main.cpp -o $@ $(LIBS) 

