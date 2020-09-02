#!/bin/bash -e

PREFIX=$CONDA_PREFIX
echo $PREFIX

# Build statically linked binary with Rust
C_INCLUDE_PATH=$PREFIX/include \
LIBRARY_PATH=$CONDA_PREFIX/lib \
LIBCLANG_PATH=$PREFIX/lib/libclang.so \
cargo build --release

# Install the binaries
C_INCLUDE_PATH=$PREFIX/include \
LIBRARY_PATH=$PREFIX/lib \
LIBCLANG_PATH=$PREFIX/lib/libclang.so \
cargo install --force --root $PREFIX

# move Rscript
cp src/snp_density_plots.R $CONDA_PREFIX/bin/