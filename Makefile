.PHONY: all help clean test clang-format

DEPS= $(wildcard src/*.h) $(wildcard src/*.cpp) $(wildcard src/bwa/.*.h) $(wildcard src/bwa/*.c)

all: ms-index

ms-index: $(DEPS)
	$(MAKE) -C src

test:
	$(MAKE) -C tests

format:
	clang-format -verbose -i src/mask.h src/index.h src/parser.h src/*.cpp tests/*.h tests/*.cpp

clean: ## Clean
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	rm -f ms-index
