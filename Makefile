.PHONY: all help clean test clang-format all-ci-linux all-ci-macos test-ci-linux test-ci-macos

DEPS= $(wildcard src/*.h) $(wildcard src/*.cpp) $(wildcard src/bwa/.*.h) $(wildcard src/bwa/*.c)

all: ms-index

ms-index: $(DEPS)
	$(MAKE) -C src

test:
	$(MAKE) -C tests

all-ci-linux:
	$(MAKE) -C src INCLUDES-PATH=/home/runner

test-ci-linux:
	$(MAKE) -C tests INCLUDES-PATH=/home/runner

all-ci-macos:
	$(MAKE) -C src INCLUDES-PATH=/Users/runner

test-ci-macos:
	$(MAKE) -C tests INCLUDES-PATH=/Users/runner

format:
	clang-format -verbose -i src/*.h src/*.cpp tests/*.h tests/*.cpp

clean: ## Clean
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	rm -f ms-index
