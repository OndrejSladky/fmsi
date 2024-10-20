#!/usr/bin/env bash
cat <(echo "#define VERSION \"") <(git describe --abbrev=4 --dirty --always --tags 2> /dev/null || cat src/version) <(echo \") | tr -d '\n' > src/version.h
