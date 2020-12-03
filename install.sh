#!/bin/bash -e

PREFIX=$CONDA_PREFIX
echo $PREFIX

# Build statically linked binary with Rust
C_INCLUDE_PATH=$CONDA_PREFIX/include \
LIBRARY_PATH=$CONDA_PREFIX/lib \
LIBCLANG_PATH=$CONDA_PREFIX/lib/libclang.so \
cargo build --release

# Install the binaries
C_INCLUDE_PATH=$PREFIX/include \
LIBRARY_PATH=$PREFIX/lib \
LIBCLANG_PATH=$PREFIX/lib/libclang.so \
cargo install --force --root $PREFIX

# Install flock
cd flock/ && pip install . && cd ../

# Install spacegraphcats
cd spacegraphcats/ && pip install . && cd ../

# move Rscript and python
# cp src/bin/snp_density_plots.R $CONDA_PREFIX/bin/
#cp src/bin/cluster.py $CONDA_PREFIX/bin/
