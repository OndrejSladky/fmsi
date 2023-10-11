.PHONY: all help clean test clang-format

SHELL=/usr/bin/env bash -eo pipefail
IND=./prophex

.SECONDARY:

DEPS= $(wildcard src/*.h) $(wildcard src/*.cpp) $(wildcard src/bwa/.*.h) $(wildcard src/bwa/*.c)

all: ms-index

ms-index: $(DEPS)
	$(MAKE) -C src

format:
	clang-format -verbose -i src/*.h src/*.cpp

clean: ## Clean
	$(MAKE) -C src clean
	rm -f ms-index
