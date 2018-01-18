# TODO: replace saraedum with sagemath eventually
ARG ARTIFACT_BASE=saraedum/sagemath-dev:latest

################################################################################
# Image containing the run-time dependencies for Sage                          #
################################################################################
FROM ubuntu:xenial as run-time-dependencies
LABEL maintainer="Erik M. Bray <erik.bray@lri.fr>"
# Set sane defaults for common environment variables.
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV SHELL /bin/bash
# Create symlinks for sage and sagemath - we copy a built sage to the target of these symlinks later.
ARG SAGE_ROOT=/home/sage/sage
RUN ln -s "$SAGE_ROOT/sage" /usr/bin/sage
RUN ln -s /usr/bin/sage /usr/bin/sagemath
# Sage needs the fortran libraries at run-time because we do not build gfortran
# with Sage but use the system's.
# We also install sudo for the sage user, see below.
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends gfortran sudo \
    && apt-get clean \
    && rm -r /var/lib/apt/lists/*
# Sage refuses to build as root, so we add a "sage" user that can sudo without a password.
# We also want this user at runtime as some commands in sage know about the user that was used during build.
ARG HOME=/home/sage
RUN adduser --quiet --shell /bin/bash --gecos "Sage user,101,," --disabled-password --home "$HOME" sage \
    && echo "sage ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-sage \
    && chmod 0440 /etc/sudoers.d/01-sage
# Run everything from now on as the sage user in sage's home
USER sage
ENV HOME $HOME
WORKDIR $HOME

################################################################################
# Image containing everything so that a make in a clone of the Sage repository #
# completes without errors                                                     #
################################################################################
FROM run-time-dependencies as build-time-dependencies
# Install the build time dependencies & git
RUN sudo apt-get update -qq \
    && sudo apt-get install -y wget build-essential automake m4 dpkg-dev sudo python libssl-dev git \
    && sudo apt-get clean \
    && sudo rm -r /var/lib/apt/lists/*
# Make everything in the build use multiple cores (this causes trouble for some packages outside Sage but it still seems to be they way Sage is doing this.)
ARG MAKE="make -j2"
ENV MAKE $MAKE

################################################################################
# Image with an empty git repository in $SAGE_ROOT.                            #
################################################################################
FROM build-time-dependencies as source-clean
ARG SAGE_ROOT=/home/sage/sage
RUN mkdir -p "$SAGE_ROOT"
WORKDIR $SAGE_ROOT
RUN git init
RUN git remote add trac git@trac.sagemath.org:sage.git
WORKDIR $HOME

################################################################################
# Image with a copy of the latest develop build in $SAGE_ROOT. We will later   #
# upgrade this copy to the contents of $HOME/context to massively speed up     #
# the build. You can disable this step to get a clean build, see below.        #
# This blows up the size of this docker image significantly, but we only use   #
# this image to create artifacts for our final image.                          #
# For a clean build, set ARTIFACT_BASE to "source-clean".                      #
################################################################################
FROM $ARTIFACT_BASE as artifact-base
FROM build-time-dependencies as source-from-latest
ARG SAGE_ROOT=/home/sage/sage
# We do not add a switch to disable this step even if we do not need it below
# as the result should usually be cached.
COPY --chown=sage:sage --from=artifact-base $SAGE_ROOT $SAGE_ROOT

################################################################################
# Image with the build context added, i.e., the directory from which `docker   #
# build` has been called in a separate directory so we can copy files from     #
# there.                                                                       #
# This blows up the size of this docker image significantly, but we only use   #
# this image to create artifacts for our final image.                          #
################################################################################
FROM source-from-latest as source-from-context
COPY --chown=sage:sage . sage-context
# Checkout the commit that is in $HOME/sage-context
ARG SAGE_ROOT=/home/sage/sage
WORKDIR $SAGE_ROOT
RUN git fetch "$HOME/sage-context" HEAD
RUN git checkout FETCH_HEAD
# Copy over all the untracked/staged/unstaged changes from sage-context
WORKDIR $HOME/sage-context
RUN if git status --porcelain | read CHANGES; then \
        git -c user.name=docker-build -c user.email=docker-build@sage.invalid stash -u \
        && git stash show -p > "$HOME"/sage-context.patch; \
    else \
        touch "$HOME"/sage-context.patch; \
    fi
WORKDIR $SAGE_ROOT
RUN patch < "$HOME"/sage-context.patch

################################################################################
# Image with a built sage but without sage's documentation.                    #
################################################################################
FROM source-from-context as make-build
# Make sure that the results runs on most CPUs.
ENV SAGE_FAT_BINARY yes
# Just to be sure Sage doesn't try to build its own GCC (even though
# it shouldn't with a recent GCC package from the system and with gfortran)
ENV SAGE_INSTALL_GCC no
RUN make build \
    && rm -r src/build/temp.*

################################################################################
# Image with a full build of sage and its documentation.                       #
################################################################################
FROM make-build as make-all
RUN make

################################################################################
# A releasable (relatively small, but still huge) image of this build with all #
# the build artifacts intact so developers can make changes and rebuild        #
#  quickly                                                                     #
################################################################################
FROM build-time-dependencies as dev
ARG SAGE_ROOT=/home/sage/sage
COPY --chown=sage:sage --from=make-all $SAGE_ROOT $SAGE_ROOT

################################################################################
# Image with a full build of sage, ready to release, i.e., with stripped       #
# binaries and some extras to run the jupyter notebook.                        #
################################################################################
FROM make-all as make-release
RUN sage -pip install terminado "notebook>=5" "ipykernel>=4.6"
RUN sage -i gap_jupyter singular_jupyter pari_jupyter

RUN export LC_ALL=C \
    && find local/lib local/bin -type f -exec strip '{}' ';' 2>&1 | grep -v "File format not recognized" |  grep -v "File truncated" || true
RUN rm -r local/share/doc/sage/doctrees local/share/doc/sage/inventory

################################################################################
# A releasable (relatively small) image which contains a copy of sage without  #
# temporary build artifacts (entrypoint: the command line interface.)          #
################################################################################
FROM run-time-dependencies as sage
ARG SAGE_ROOT=/home/sage/sage
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/sage $SAGE_ROOT/sage
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/local $SAGE_ROOT/local
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/src/bin $SAGE_ROOT/src/bin
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/VERSION.txt $SAGE_ROOT/VERSION.txt
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/build $SAGE_ROOT/build
# Put scripts to start gap, gp, maxima, ... in /usr/bin
RUN sudo sage --nodotsage -c "install_scripts('/usr/bin')"
COPY ./docker/entrypoint/sage-cli.sh /usr/local/bin/sage-entrypoint
ENTRYPOINT /usr/local/bin/sage-entrypoint

################################################################################
# A releasable (relatively small) image which contains a copy of sage without  #
# temporary build artifacts (entrypoint: the jupyter notebook.)                #
################################################################################
FROM sage as sage-jupyter
EXPOSE 8888
COPY ./docker/entrypoint/jupyter.sh /usr/local/bin/sage-entrypoint
