#!/bin/bash -e

PREFIX=$CONDA_PREFIX
echo $PREFIX
VERSION=0.6.0

# copy statically linked binary
cp releases/v$VERSION/* $CONDA_PREFIX/bin/

# Install flight
cd flight/ && pip install . && cd ../