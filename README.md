# Workflow: Estimate VerifyBamID2 Reference Data

This Snakemake workflow estimates parameter files for a new genome build for [VerifyBamID2](https://github.com/Griffan/VerifyBamID).
It requires a reference genome fasta and one or more variant call files, anticipated to be from the 1000 Genomes Project.

## Authors

* Lightning Auriga (@lightning-auriga)

## Usage

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@github.com:lightning-auriga/estimate-verifybamid2-reference-data.git
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
in order of chromosome number and further sorted by physical position within chromosome. This is a restriction of
[bcftools concat](https://samtools.github.io/bcftools/bcftools.html#concat), and while annoying, it is easily obtained from
existing 1000 Genomes callsets by simply organizing the vcfs within the manifest in chromosome order. If this poses some sort
of problem for any future issue, please go ahead and open an issue, and additional sort operations can be patched in.

For cluster submission of jobs, additional tool-specific resource specifications are exposed in `config/config_resources.yaml`.
The default settings for this workflow are configured based on a Slurm scheduler on an virtual on-demand HPC using
[AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html). At the time of writing,
all tasks except for the parameter estimation with VerifyBamID2 itself can run on a single threaded t3.medium node with 4G RAM;
run on a spot fleet, these tasks are very affordable. The Verify run itself requires additional RAM, and is currently configured to run
on a c7i.8xlarge node; a smaller node would certainly be functional as well.

The most intensive aspect of this workflow is really disk space, when downloading the full vcfs from e.g. the 1000 Genomes Project.
The workflow is confirmed to run on an FSx for Lustre drive with 1.2T of available space (the minimum possible request for such a partition).
The combined cost for this workflow on the above-configured EC2 resources is approximately $5.00 USD.

Note that resource configuration is specified for a slurm job scheduler using a slurm cluster profile. The resources are specified
as the `slurm_partition` parameter in per-rule resource configuration. If you are interested in running this workflow with a profile
for a different scheduler, you can either adapt the profile to recognize the `slurm_partition` parameter, or edit `slurm_partition` to
whatever parameter (e.g. `queue`) your profile is expecting. For the latter case, the command `sed -i 's/slurm_partition/queue/' workflow/rules/*smk`
will do the trick.

### Step 3: Install Snakemake

Install Snakemake using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html):

    mamba create -c conda-forge -c bioconda -n snakemake snakemake

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

The output VerifyBamID2 SVD estimate files are emitted to `results/verify_output` with the prefix specified
by the user in `config["dataset-name"]`. The relevant file extensions are `*.UD`, `*.V`, `*.mu`, and `*.bed`;
see [the VerifyBamID2 GitHub page](https://github.com/Griffan/VerifyBamID) for more documentation.

After the pipeline runs to completion, executing `snakemake -j1 --use-conda -p benchmark_report` will create an html
in `results/performance_benchmarks/performance_benchmarks.html`, describing per-job runtime performance.
The statistics reported are unfamiliar to many;
see [this remarkably helpful thread](https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake)
for more documentation if interested. This report is only interesting if the user wishes to try to optimize task performance.
