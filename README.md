**Institut Curie - Nextflow scChIPseq-inDrop analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow manager to run tasks across multiple compute infrastructures in a very portable manner.
It supports [conda](https://docs.conda.io) package manager and  [singularity](https://sylabs.io/guides/3.6/user-guide/) / [Docker](https://www.docker.com/) containers making installation easier and results highly reproducible.

## Pipeline summary

This data engineering pipeline is designed to treat single-cell chromatin Immuno-Precipitation sequencing from raw reads (fastq, paired end) to exploitable count matrix. The multiple steps involved in the pipeline are :

1. Align R2 reads on barcode indexes
2. Trimming R2 reads
3. Align paired-end reads on the genome
4. Assignation of cell barcodes to aligned read
5. Removal of Reverse Transcription (RT) & Polymerase Chain Reaction (PCR) duplicates
6. Removal of reads based on window screening (if Read2 was unmapped)
7. Counting (Generation of count matrix)
8. Generation of coverage file (bigwig) (CPM normalization)
9. Reporting

![MultiQC](docs/images/Simplified_pipeline.png =x150)


### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
=======================================================

Usage:

Mandatory arguments:
--reads [file]                   Path to input data (must be surrounded with quotes)
--samplePlan [file]              Path to sample plan file if '--reads' is not specified
--genome [str]                   Name of the reference genome. See the `--genomeAnnotationPath` to defined the annotation path
-profile [str]                   Configuration profile to use (multiple profiles can be specified with comma separated values)

Skip options: All are false by default
--skipSoftVersion [bool]         Do not report software version
--skipMultiQC [bool]             Skip MultiQC

Other options:
--outDir [dir]                  The output directory where the results will be saved
-w/--work-dir [dir]             The temporary directory where intermediate data will be saved
-name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
-oldDesign [bool]               If dark cycles are used write this option as true. Default is false.
-keepRTdup [bool]               Keep RT duplicats. Default is false.
-window [int]                   Select the window size. Default is 50.
-minCounts [int]                Select the minimum count per barcodes after removing duplicates. Default is 1000.
-removeBlackRegion [bool]       Remove black region. Default is true.
-mark [str]                     Histone mark targeted 'h3k27me3', 'h3k4me3' or 'unbound'. Default is 'h3k27me3'.
-binSize [int]                  Bin size to use (in base pairs). Default is 50000. 

=======================================================
Available profiles
-profile test                    Run the test dataset
-profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile path                    Use the installation path defined for all tools. Use `--globalPath` to define the installation path
-profile multipath               Use the installation paths defined for each tool. Use `--globalPath` to define the installation path
-profile docker                  Use the Docker images for each process
-profile singularity             Use the Singularity images for each process. Use `--singularityPath` to define the path of the singularity containers
-profile cluster                 Run the workflow on the cluster, instead of locally

```


### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follows:

#### Run the pipeline on a test dataset

See the file `conf/test.config` to set your test dataset.

```bash
nextflow run main.nf -profile test,conda

```

#### Run the pipeline from a `sample plan` and a `design` file

```bash
nextflow run main.nf --samplePlan mySamplePlan.csv --design myDesign.csv --genome 'hg19' --genomeAnnotationPath /my/annotation/path --outDir /my/output/dir

```

### Defining the '-profile'

By default (whithout any profile), Nextflow excutes the pipeline locally, expecting that all tools are available from your `PATH` environment variable.

In addition, several Nextflow profiles are available that allow:
* the use of [conda](https://docs.conda.io) or containers instead of a local installation,
* the submission of the pipeline on a cluster instead of on a local architecture.

The description of each profile is available on the help message (see above).

Here are a few examples to set the profile options:

#### Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
```bash
-profile path --globalPath /my/path/to/bioinformatics/tools
```

#### Run the pipeline on the cluster, using the Singularity containers
```bash
-profile cluster,singularity --singularityPath /my/path/to/singularity/containers
```

#### Run the pipeline on the cluster, building a new conda environment
```bash
-profile cluster,conda --condaCacheDir /my/path/to/condaCacheDir

```

For details about the different profiles available, see [Profiles](docs/profiles.md).

### Sample plan

A sample plan is a csv file (comma separated) that lists all the samples with a biological IDs.
The sample plan is expected to contain the following fields (with no header):

```
SAMPLE_ID,SAMPLE_NAME,path/to/R1/fastq/file,path/to/R2/fastq/file (for paired-end only)
```

### Genome annotations

The pipeline does not provide any genomic annotations but expects them to be already available on your system. The path to the genomic annotations can be set with the `--genomeAnnotationPath` option as follows:

```bash
nextflow run main.nf --samplePlan mySamplePlan.csv --design myDesign.csv --genome 'hg19' --genomeAnnotationPath /my/annotation/path --outDir /my/output/dir

```

For more details see  [Reference genomes](docs/referenceGenomes.md).

## Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/referenceGenomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Credits

This pipeline has been written by the DEPIC team, the single cell platform & the bioinformatics platform of the Institut Curie (Pacôme Prompsy, Louisa Hadj Abed, Celine Vallot, Nicolas Servant)

## Contacts

For any question, bug or suggestion, please use the issue system or contact the bioinformatics core facility.
