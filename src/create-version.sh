#!/usr/bin/env bash

# Write version.h: the version of the file `version`, followed by the commit of the working tree if there is one
# (and "-dirty" if it has changes that are not committed), e.g. v2.0.0-dev-gbc4cc31. The version does not depend
# on tags, so that a development version is not named after the last release.

cd "$(dirname "$0")" || exit 1

version="v$(tr -d '[:space:]' < version)"
commit=$(git rev-parse --short=7 HEAD 2> /dev/null)
if [ -n "$commit" ]; then
    version="$version-g$commit"
    if ! git diff --quiet HEAD -- . ../tests 2> /dev/null; then
        version="$version-dirty"
    fi
fi
printf '#define VERSION "%s"' "$version" > version.h
