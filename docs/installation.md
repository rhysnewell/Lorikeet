---
title: Installation
---

Installation
========

#### Option 1: Build manually
You may need to manually set the paths for `C_INCLUDE_PATH`, `LIBRARY_PATH`, `LIBCLANG_PATH`, and `OPENSSL_DIR` to their corresponding
paths in the your conda environment if they can't properly be found on your system. This method also assumes you have 
previously installed rust via rustup on your system.
```
GIT_LFS_SKIP_SMUDGE=1 git clone --recursive https://github.com/rhysnewell/Lorikeet.git;
cd Lorikeet;
conda env create -n lorikeet -f lorikeet.yml; 
conda activate lorikeet;
pip install --upgrade cmake;
bash install.sh # or run without installing e.g. `cargo run --release -- call -h`;
lorikeet call -h
```

Depending on your local network configuration, you may have problems obtaining Lorikeet via git.
If you see something like this you may be behind a proxy that blocks access to standard git:// port (9418).

```
$ git clone --recursive git://github.com/rhysnewell/Lorikeet.git
Cloning into 'Lorikeet'...
fatal: Unable to look up github.com (port 9418) (Name or service not known)
```

Luckily, thanks to this handy tip from the developer of [Freebayes](https://github.com/ekg/freebayes) we can work around it.
If you have access to https:// on port 443, then you can use this 'magic' command as a workaround to enable download of the submodules:

```
git config --global url.https://github.com/.insteadOf git://github.com/
```

#### Option 2: Conda 

Install into current conda environment:
```
conda install -c conda-forge bioconda lorikeet-genome
```

Create fresh conda environment and install lorikeet there:
```
conda create -n lorikeet -c conda-forge bioconda lorikeet-genome && \
conda activate lorikeet
```

#### Option 3: Install static binary **Not currently recommended**
The static binary is the easiest to use, however it is compiled with `musl` which for some reason makes Lorikeet (and other
rust binaries) perform much slower than they usually should. As such, we do not recommend using this static binary until
this problem is sorted out.

You can make use of the precompiled static binaries that come with this repository. You will have to install the lorikeet
conda environment using the lorikeet.yml.
```
GIT_LFS_SKIP_SMUDGE=1 git clone --recursive https://github.com/rhysnewell/Lorikeet.git;
cd Lorikeet;
conda env create -n lorikeet -f lorikeet.yml;
conda activate lorikeet
```

Once you have created the conda environment download and install the latest release file from github
```
wget https://github.com/rhysnewell/Lorikeet/releases/download/latest/lorikeet-x86_64-unknown-linux-musl-v0.6.1.tar.gz;
tar -xvzf lorikeet-x86_64-unknown-linux-musl-v*.tar.gz;
cp release/lorikeet $CONDA_PREFIX/bin;
cp release/remove_minimap2_duplicated_headers $CONDA_PREFIX/bin;
```

