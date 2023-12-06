---
title: Installation
---

Installation
========

### Installation

Lorikeet is distributed via Crates.io https://crates.io/crates/lorikeet-genome. Additional packages can be downloaded via conda using the `lorikeet.yml` file provided. Ensure that cargo is installed on your system:

```bash
curl https://sh.rustup.rs -sSf | sh
```

Then install lorikeet:

```bash
cargo install lorikeet-genome
```

Alongside required packages:

```bash
conda env create -f lorikeet.yml -n lorikeet
conda activate lorikeet
```