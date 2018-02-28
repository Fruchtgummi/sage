#!/bin/sh

# This script gets called from CI to push a docker image for binder to Docker
# Hub and also to push corresponding binder configuration to the github
# repository $SAGE_BINDER_ENV_GITHUB.

# The following environment variables are required in your CI system:
# * DOCKER_USER (your username on Docker Hub)
# * SECRET_DOCKER_PASS (your password on Docker Hub)
# * SAGE_GITHUB (your Sage repository on github, e.g., username/sagemath)
# * SAGE_BINDER_ENV_GITHUB (a repository on github for autogenerated binder
#   configurations, e.g., username/sage-binder-env)
# * SECRET_SAGE_BINDER_ENV_KEY (an RSA private key with push access to
#   SAGE_BINDER_ENV_GITHUB)

# Warning: It is a bit scary to give CI push access to github. We try to make
# sure that your key does not leak into the logs, see ./protect-secrets.sh.
# However, please make sure that the key has only push access to
# SAGE_BINDER_ENV_GITHUB and that you do nothing otherwise too important in
# that repository.

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

if [[ -z "$DOCKER_TAG" || -z "$DOCKER_USER" || -z "$SECRET_DOCKER_PASS" ]]; then
    echo "Build can not be pushed to docker hub. Not pushing binder configuration."
    exit 0
fi

if [[ -z "$SECRET_SAGE_BINDER_ENV_KEY" ]]; then
    echo "No deployment key for github configured. Not pushing binder configuration."
    exit 0
fi

if [[ -z "$SAGE_BINDER_ENV_GITHUB" ]]; then
    echo "No github repository configured for binder builds. Not pushing binder configuration."
    exit 0
fi

# Tag the sagemath:$DOCKER_TAG image that CI has just built as
# sagemat-review:COMMIT_HASH so we can refer to it uniquely later.
[[ "$DOCKER_TAG" = "master" ]] && DOCKER_TAG=latest
export COMMIT_HASH=$(git log -1 --format=format:%h)
export DOCKER_IMAGE="${DOCKER_USER}/sagemath-review:${COMMIT_HASH}"
docker tag "${DOCKER_USER}/sagemath:${DOCKER_TAG}" "$DOCKER_IMAGE"
cat "$SECRET_DOCKER_PASS" | docker login -u $DOCKER_USER --password-stdin
docker push "$DOCKER_IMAGE"

# Restore private key for ${SAGE_BINDER_ENV_GITHUB}.
mkdir -p ~/.ssh/
cat "$SECRET_SAGE_BINDER_ENV_KEY" | sed 's/\\n/\n/g' > ~/.ssh/id_rsa
chmod 700 ~/.ssh
chmod 600 ~/.ssh/id_rsa

# Collect some metadata to include in the home page of the Jupyter notebook and
# also in the README of the branch on SAGE_BINDER_ENV_GITHUB.

# json_escape() from https://stackoverflow.com/a/13466143/812379
json_escape() {
    printf '%s' $1 | python -c 'import json,sys; print(json.dumps(sys.stdin.read()))'
}
export AUTHOR=$(json_escape "`git log -1 --format=format:%an`")
export COMMIT_MESSAGE=$(json_escape "`git log -1 --format=format:%s%n%n%-b`")
export COMMIT_TIMESTAMP=$(json_escape "`git log -1 --format=format:%aD`")
export COMMIT_URL="https://github.com/${SAGE_GITHUB:-sagemath/sage}/commit/$(git log -1 --format=%H)"
export BINDER_URL="https://mybinder.org/v2/gh/${SAGE_BINDER_ENV_GITHUB}/${BRANCH}?filepath=review.ipynb"

# Substitute the above variables in all files in .ci/binder.
cd .ci/binder
git init
for template in *;do
    mv "$template" "${template}.tmpl"
    envsubst < "${template}.tmpl" > "$template"
    git add "$template"
done
# Verify that the notebook is valid JSON
python -m json.tool < review.ipynb > /dev/null

# Force push a new README and configuration to BRANCH on SAGE_BINDER_ENV_GITHUB.
git -c user.name=circleci -c user.email=circleci@build.invalid commit -m "automatically generated from template"
unset SSH_AUTH_SOCK
unset SSH_ASKPASS
git push --force "git@github.com:${SAGE_BINDER_ENV_GITHUB}.git" "HEAD:${BRANCH}"
cd ../..
