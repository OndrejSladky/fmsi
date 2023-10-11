CC=			g++
#CC=			clang --analyze
CFLAGS=		-g -Wall -Wno-unused-function -O2 -g3
RAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)

PROG=../ms-index

AOBJS2=	\
			bwa/bwashm.o \
			bwa/bwase.o \
			bwa/bwaseqio.o \
			bwa/bwtgap.o \
			bwa/bwtaln.o \
			bwa/bamlite.o \
			bwa/bwape.o \
			bwa/kopen.o \
			bwa/pemerge.o \
			bwa/maxk.o \
			bwa/bwtsw2_core.o \
			bwa/bwtsw2_main.o \
			bwa/bwtsw2_aux.o \
			bwa/bwt_lite.o \
			bwa/bwtsw2_chain.o \
			bwa/fastmap.o \
			bwa/bwtsw2_pair.o \
			bwa/QSufSort.o \


INCLUDES=	-Ibwa
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.cpp.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


all:$(PROG)

clean:
	rm -f gmon.out *.o a.out $(PROG) *~ *.a
	# if BWA Makefile is present
	test -f bwa/Makefile && $(MAKE) -C bwa clean

$(PROG): bwa/libbwa.a $(AOBJS2) bwa/QSufSort.c $(wildcard *.cpp *.h)
	$(CC) $(INCLUDES) $(CFLAGS) $(DFLAGS) $(AOBJS2) bwa/QSufSort.c main.cpp -o $@ -Lbwa -lbwa $(LIBS)

#bwa/libbwa.a $(AOBJS2) bwtexk.o:
bwa/libbwa.a:
	 $(MAKE) -C bwa
