#!/usr/bin/env bash
cat <(echo "#define VERSION \"") <(git describe --abbrev=4 --dirty --always --tags 2> /dev/null || cat version) <(echo \") | tr -d '\n' > version.h
