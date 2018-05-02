#!/bin/sh

# This script gets called from CI to push a docker image for binder to Docker
# Hub and also to push corresponding binder configuration to the github
# repository $SAGE_BINDER_ENV_GITHUB.

# The following environment variables are required in your CI system:
# * DOCKER_IMAGE the public name of the sagemath image for which we should
#   build a binder image
# * HTTP_GIT_SAGE (used to build a link to the commit represented by this
#   binder image as $HTTP_GIT_SAGE/COMMIT_SHA )
# * SSH_GIT_BINDER (repository to push the autogenerated binder configuration
#   to, git@ URL)
# * HTTP_GIT_BINDER (repository to pull the autogenerated binder configuration
#   from, https:// URL)
# * SECRET_SSH_GIT_BINDER_KEY (an RSA private key with push access to GIT_BINDER)

# Warning: It is a bit scary to give CI push access to git. We try to make sure
# that your key does not leak into the logs, see ./protect-secrets.sh.
# However, please make sure that the key has only push access to GIT_BINDER and
# that you do nothing otherwise too important in that repository.

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

if [[ -z "$DOCKER_IMAGE" ]]; then
    echo "DOCKER_IMAGE is not set. The image is not available on a public registry. Cannot build binder."
    exit 0
fi

if [[ -z "$SECRET_GIT_BINDER_KEY" ]]; then
    echo "No deployment key for git configured. Not pushing binder configuration."
    exit 0
fi

if [[ -z "$SSH_GIT_BINDER" ]]; then
    echo "No git repository configured for binder builds. Not pushing binder configuration."
    exit 0
fi

# Restore private key for ${SECRET_SSH_GIT_BINDER_KEY}.
mkdir -p ~/.ssh/
cat "$SECRET_SSH_GIT_BINDER_KEY" | sed 's/\\n/\n/g' > ~/.ssh/id_rsa
chmod 700 ~/.ssh
chmod 600 ~/.ssh/id_rsa

# escape_md: Escape for interpolation in Markdown literal blocks by stripping out all backticks.
escape_md() {
    echo -n "$1" | sed 's/`//g'
}
# escape_json: Escape for interpolation in JSON double quoted strings.
escape_json() {
    echo -nE "$1" | python -c 'import json,sys; print(json.dumps(sys.stdin.read())[2:-2])'
}
# Collect some metadata to include in the home page of the Jupyter notebook and
# also in the README of the branch on SSH_GIT_BINDER.
export AUTHOR_MD=$(escape_md "`git log -1 --format=format:%an`")
export AUTHOR_JSON=$(escape_json "$AUTHOR_MD")
export COMMIT_MESSAGE_MD=$(escape_md "`git log -1 --format=format:%s%n%n%-b`")
export COMMIT_MESSAGE_JSON=$(escape_json "$COMMIT_MESSAGE_MD")
export COMMIT_TIMESTAMP=$(git log -1 --format=format:%aD)
export COMMIT_URL="${HTTP_GIT_SAGE}/$(git log -1 --format=%H)"
export BINDER_URL="https://mybinder.org/v2/git/${HTTP_GIT_BINDER}/${BRANCH}?filepath=review.ipynb"

# Substitute the above variables in all files in .ci/binder.
cd .ci/binder
git init
for template in *;do
    mv "$template" "${template}.tmpl"
    envsubst < "${template}.tmpl" > "$template"
    rm -f "${template}.tmpl"
    git add "$template"
done
# Verify that the notebook is valid JSON
python -m json.tool < review.ipynb > /dev/null

# Force push a new README and configuration to BRANCH on SSH_GIT_BINDER.
git -c user.name=sage.binder -c user.email=sage.binder@build.invalid commit -m "automatically generated from template"
unset SSH_AUTH_SOCK
unset SSH_ASKPASS
git push --force "${SSH_GIT_BINDER}" "HEAD:${BRANCH}"
