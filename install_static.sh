#!/bin/bash -e

PREFIX=$CONDA_PREFIX
echo $PREFIX

# copy statically linked binary
cp releases/* $CONDA_PREFIX/bin/

# Install flight
cd flight/ && pip install . && cd ../