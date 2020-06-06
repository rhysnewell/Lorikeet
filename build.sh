#!/bin/bash -e

echo $CONDA_PREFIX

# Build statically linked binary with Rust
C_INCLUDE_PATH=$CONDA_PREFIX/include \
LIBRARY_PATH=$CONDA_PREFIX/lib \
LIBCLANG_PATH=$CONDA_PREFIX/lib/libclang.so \
cargo build --release

# Install the binaries
cargo install --root $CONDA_PREFIX