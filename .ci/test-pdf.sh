#!/bin/sh

# This script gets called from CI to check that the PDF documentation builds without errors.

# Usage: ./test-pdf.sh sage-dev-image

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -exo pipefail

cat << EOF | docker run -i "$1" bash || exit 1
set -ex
sudo apt-get update
sudo apt-get install -y --no-install-recommends latexmk texlive-latex-recommended texlive-fonts-recommended texlive-latex-extra texlive-generic-extra
sage -docbuild all pdf
EOF
