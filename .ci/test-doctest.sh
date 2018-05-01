#!/bin/sh

# This script gets called from CI to run doctests in the sagemath build

# Usage: ./test-doctest.sh IMAGE-NAME --new|--short|--long

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -ex

run_doctests() {
    docker run --entrypoint sh "$1" -c "sage -tp $@ ||
                                        sage -tp --failed $@ ||
                                        sage -tp --failed $@"
}

case "$2" in
    --new)
        LATEST_RELEASE=`git log --author release@sagemath.org -1 --format=%H`
        MODIFIED_ADDED=`git diff --diff-filter=MC --name-only "$LATEST_RELEASE" | grep -E '*.(py|pyx|rst)'`
        if [[ "x$MODIFIED_ADDED" = "x" ]]; then
            echo "No testable files have been modified/added."
        else
            run_doctests $MODIFIED_ADDED
        fi
        ;;
    --short)
        # TODO: Upgrade this with https://trac.sagemath.org/ticket/25270
        run_doctests --all
        ;;
    --long)
        run_doctests --all
        ;;
    *)
        exit 1
        ;;
esac
