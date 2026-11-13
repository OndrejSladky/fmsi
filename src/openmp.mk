# OpenMP is used by libsais only, which then constructs the suffix array on several threads
# (fmsi index2 -t); the rest of FMSI is threaded with std::thread and does not need it.
#
# OpenMP is switched on if the C compiler supports it:
#  - GCC and LLVM clang take -fopenmp (Linux, conda, and Homebrew GCC on macOS:
#    make CC=gcc-16 CXX=g++-16);
#  - Apple clang has no OpenMP runtime of its own and takes the one of Homebrew
#    (brew install libomp), found through LIBOMP.
# Otherwise, or with make OPENMP=0, the suffix array is constructed on one thread.
# CC and CXX must be compilers of the same family, as they must agree on the OpenMP runtime.

OPENMP ?= 1
OPENMP_CFLAGS :=
OPENMP_LIBS :=

ifneq ($(OPENMP),0)
  OPENMP_PROBE := $(dir $(lastword $(MAKEFILE_LIST)))openmp_probe.c
  # openmp_works <flags>: yes if the probe compiles with the flags.
  openmp_works = $(shell $(CC) $(1) -fsyntax-only $(OPENMP_PROBE) 2> /dev/null && echo yes)
  ifeq ($(call openmp_works,-fopenmp),yes)
    OPENMP_CFLAGS := -fopenmp -DLIBSAIS_OPENMP
    OPENMP_LIBS := -fopenmp
  else
    LIBOMP ?= $(shell brew --prefix libomp 2> /dev/null)
    ifneq ($(wildcard $(LIBOMP)/lib/libomp.*),)
      ifeq ($(call openmp_works,-Xclang -fopenmp -I$(LIBOMP)/include),yes)
        OPENMP_CFLAGS := -Xclang -fopenmp -I$(LIBOMP)/include -DLIBSAIS_OPENMP
        OPENMP_LIBS := -L$(LIBOMP)/lib -lomp
      endif
    endif
  endif
endif

# The setting is kept in a file that the objects depend on, so that changing it rebuilds them.
OPENMP_SETTING := $(dir $(lastword $(MAKEFILE_LIST)))openmp.setting
ifneq ($(strip $(shell cat $(OPENMP_SETTING) 2> /dev/null)),$(strip $(CC) $(OPENMP_CFLAGS)))
  $(shell echo '$(strip $(CC) $(OPENMP_CFLAGS))' > $(OPENMP_SETTING))
endif

ifeq ($(OPENMP_CFLAGS),)
  $(info OpenMP: off (the suffix array is constructed on one thread))
else
  $(info OpenMP: on ($(OPENMP_CFLAGS)))
endif
