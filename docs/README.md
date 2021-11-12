![](/images/lorikeet_logo.png)

![](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/license.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/version.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/platforms.svg)


Lorikeet is a within-species variant analysis pipeline for metagenomic communities that utilizes both long and short reads.
Lorikeet utilizes a re-implementaion of the GATK HaplotypeCaller algorithm, performing local re-assembly of potentially active
regions within candidate genomes. Called variants can be clustered into likely strains using a combination of UMAP and HDBSCAN.
Additional statistics, like consensus ANI, population ANI, and subpopulation ANI will also be calculated for each input
geome providing values for each sample compared to the reference and also compared to all other samples.

Lorikeet has a variety of subcommands with the main being `call` and `genotype`. The `call` pipeline will take any number
of input genomes and samples and perform robust variant calling and ANI calculations. The `genotype` algorithm takes this
a step further and attempts to reconstruct strain haplotypes from the called variants and return complete strain genomes.

## Let's add another page

All of your documentation lives under the `docs` directory. You can start adding markdown files, and
when running the `serve` command you will see changes automatically updated in the browser.

Try it - add a new markdown file under `docs`, paste the content below, and watch what happens.

```markdown
---
title: Another page
---

Adding new pages is that simple
===============================

```

When you hit save, you should see the left side navigation has updated, and a link to your new page
shows up.

## What next?

There are plenty of resources to learn more about Doctave and how to use it effectively. Here are
some articles to get you started:

* [The official tutorial](https://cli.doctave.com/tutorial)
* [Deployment instructions](https://cli.doctave.com/deployment)
* [Doctave docs](https://cli.doctave.com/)

## Where can I get help?

Feel free to open issues on the [Github Repo](https://github.com/Doctave/doctave), especially if
you did not find an answer to a question in our documentation. You can also reach out directly to
the maintainer via [Twitter](https://twitter.com/NiklasBegley) or email at nik@doctave.com.
