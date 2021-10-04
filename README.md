[![Snakemake](https://img.shields.io/badge/snakemake-≥4.6.0-brightgreen.svg)](https://snakemake.bitbucket.io)

# U251 glioblastoma cell line data analysis pipeline

Reproducible pipeline of the experiments described in the manuscript "..."


## Requirements

* Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).
* working `heinz` executable in your path. Install from [https://github.com/ls-cwi/heinz].
* To visualize the results: `examine.jar`. Build from [https://github.com/AlBi-HHU/eXamine-stand-alone].


## Usage

### Step 1: Setup system

#### Variant a: Installing Miniconda on your system

If you are on a Linux system with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6), you can simply install Miniconda3 with

    curl -o /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/miniconda.sh

Make sure to answer `yes` to the question whether your PATH variable shall be modified.
Afterwards, open a new shell/terminal.

#### Variant b: Use a Docker container

Otherwise, e.g., on MacOS or if you don't want to modify your system setup, install [Docker](https://www.docker.com/), run

    docker run -it continuumio/miniconda3 /bin/bash

and execute all the following steps within that container.

#### Variant c: Use an existing Miniconda installation

If you want to use an existing Miniconda installation, please be aware that this is only possible if it uses Python 3 by default. You can check this via

    python --version

Further, ensure it is up to date with

    conda update --all

### Step 2: Setup Bioconda channel

Setup Bioconda with

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Step 3: Install bioconda-utils and Snakemake

Install bioconda-utils and Snakemake >=4.6.0 with

    conda install bioconda-utils snakemake

If you already have an older version of Snakemake, please make sure it is updated to >=4.6.0.

### Step 4: Download the workflow

First, create a working directory:

    mkdir u251-ivgpcr-workflow
    cd u251-ivgpcr-workflow

Then, download the workflow archive from TODO and unpack it with

    tar -xf u251-ivgpcr-workflow.tar.gz

### Step 5: Run the workflow

Execute the analysis workflow with Snakemake

    snakemake --use-conda

Please wait a few minutes for the analysis to finish.

If you have been running the workflow in the docker container (see above),
you can obtain the results with

    docker cp <container-id>:/bioconda-workflow/figs .

whith `<container-id>` being the ID of the container.


## Known errors

* If you see an error like
  ```
  ImportError: No module named 'appdirs'
  ```
  when starting Snakemake, you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with

      conda update --all

  and then reinstall the `appdirs` and `snakemake` package with

      conda install -f appdirs snakemake
* If you see an error like
  ```
  ImportError: Missing required dependencies ['numpy']
  ```
  you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with

      conda update --all

  and then reinstall the `snakemake` package with

      conda install -f snakemake

## Remarks regarding reproducibility

GO terms and KEGG pathways are snapshots from March 2018. They can be replaced by new, udated versions by deleting the corresponding files in the `GO/` and `KEGG/` subdirectories. Note, however, that the results of the enrichment process may change.
<<<<<<< HEAD
=======

same holds vor heinz.

## TODO

* Write this README
* networks --> provide as data or get via Snakemake?
  * up to date?
* start w/raw data?
* KEGG pathway enrichment -- discuss background. rule needs condas scipy, pandas
* multiple modules
* root at STAT3
* differences between modules? Venn diagram or cooler?
* How to distribute eXamine so that it runs?
* get multiple testing correction in GO and KEGG enrichment
* Rechte? Zum Beispiel an dem KEGG-Kram, den wir ins Repository gelegt haben?
* double-check enrichment. values ar sooo low for big modules. realistic?
>>>>>>> 5f440f8643e386d90a798325ad73f773daf1f28b
