.PHONY: all clean all-ci-linux all-ci-macos
CC=										g++
CFLAGS=								-g -Wall -Wno-unused-function -g3
PROG=									../ms-index
INCLUDES=							-I~/include -L~/lib
INCLUDES-CI-LINUX=		-I/home/runner/include -L/home/runner/lib
INCLUDES-CI-MACOS=		-I/Users/runner/include -L/Users/runner/lib
LIBS=									-lz -lsdsl -ldivsufsort -ldivsufsort64

all:$(PROG)

all-ci-linux:
	$(CC) $(INCLUDES-CI-LINUX) $(CFLAGS)  QSufSort.c main.cpp -o $@ $(LIBS) 
	
all-ci-macos:
	$(CC) $(INCLUDES-CI-MACOS) $(CFLAGS)  QSufSort.c main.cpp -o $@ $(LIBS) 

clean:
	rm -f $(PROG)	

$(PROG): $(wildcard *.cpp *.c *.h)
	$(CC) $(INCLUDES) $(CFLAGS) QSufSort.c main.cpp -o $@ $(LIBS) 

