---
title: 'Install Lorikeet'
date: 2019-02-11T19:27:37+10:00
weight: 2
---
### Dependencies

For read mapping:
 * [Minimap2]() by default
 * [BWA]()
 * [SAMtools]()

For open reading frame calling:
 * [Prodigal]()

Additionally, you'll need Python >= 3.6 with the following packages:
 * [Nimfa]()
 * [Numpy]()
 * [SciPy]()
 * [Threadpoolctl]()
### Easy install

A static binary of the latest release of lorikeet can be downloaded directly from [GitHub](https://github.com/rhysnewell/Lorikeet/releases/latest)

```commandline
tar -xvzf lorikeet-x86_64-unknown-linux-musl.tar.gz
```

You'll then need to place both the Lorikeet binary and the nmf.py script into your current bin/ folder

We are working on integrating the nmf.py script into Rust, however current methods using Opy3 require unstable versions of Rust.

### Build from scratch

Building from scratch requires Rust to be installed on your machine. Installation instructions for rust can be found 
[here](https://www.rust-lang.org/tools/install).

Once installed, clone the latest Lorikeet version from GitHub.

```
git clone https://github.com/rhysnewell/Lorikeet.git && cd Lorikeet/
```

