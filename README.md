# Workflow: Estimate VerifyBamID2 Reference Data

This Snakemake workflow estimates parameter files for a new genome build for [VerifyBamID2](https://github.com/Griffan/VerifyBamID).
It requires a reference genome fasta and one or more variant call files, anticipated to be from the 1000 Genomes Project.

## Authors

* Lightning Auriga (@lightning-auriga)

## Usage

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@github.com:lightning-auriga/estimate-verifybamid2-reference-data
```

### Step 2: Configure workflow

The following configuration settings can be adjusted in `config/config.yaml`:

|Setting|Description|
|---|---|
|`manifest`|path to workflow manifest; should probably remain the default `config/manifest.tsv`|
|`reference-fasta`|path to the reference genome fasta file. can be a local path, an https:// URL, or an s3:// path|
|`downsampled-variant-count`|how many high frequency, high call rate variants to downsample from the input vcfs. according to VerifyBamID2 documentation, should be at least 5000|
|`dataset-name`|name to prefix output files with|

In the manifest, by default at `config/manifest.tsv`, provide paths to vcf files containing variant calls for this reference genome.
As with the reference fasta, these can be local paths, https:// URLs, or s3:// paths. For the moment, these are expected to be
one vcf per chromosome for a combined set of 1000 Genomes Project autosome calls. Note that the files should currently be sorted
in order of chromosome number; this restriction will be relaxed when I'm less lazy.

### Step 3: Install Snakemake

Install Snakemake using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html):

    mamba create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    mamba activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --profile /path/to/cluster/profile --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

The results files will, temporarily, be in the directory `results/combine_vcfs/`, due to silly restrictions in naming
from VerifyBamID2. I'll update these to be less ridiculous moving forward. Your guess is as good as mine what the
files `*.UD`, `*.mu`, and `*.bed` contain.
