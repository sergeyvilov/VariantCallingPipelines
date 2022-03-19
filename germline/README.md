## Call germline variants with GATK

This is a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for calling germline variants according to [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).

The pipeline generates VCF files with germline variants for given BAM samples.

Input samples should be prepared according to [GATK best practices of data-preprocessing for variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)

## Setup

1. Install  [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Create and activate the conda environment

```
conda env create -f environment.yml
conda activate vcalling

```

3. Install GATK4 tools

There are two possible ways:

-  Use gatk4 package from the *vcalling* environment
In this case, not all tools can run properly (e.g. VariantRecalibrator can't make R plots)
But it's ok for such tools as Mutect2 or Haplotype Caller

- Use the native gatk package from the Broad institute
All the tools should work properly in this case.

Assuming the current folder contains the snakefile:

- Download the latest gatk release as a zip archive
https://github.com/broadinstitute/gatk/releases
- Extract the archive to some gatk path
- Add gatk execution path to the config.yaml:
  gatk_path: some_path/gatk-4.2.3.0/gatk
- cp some_path/gatk-4.2.3./gatkcondaenv.yml ./envs
- With envs/gatkcondaenv.yml:
  * Add the line
    "- openjdk=8.0"
  * In the last line,
    change gatkPythonPackageArchive.zip to some_path/gatk-4.2.3./gatkPythonPackageArchive.zip
- For all rules requiring gatk add:
  conda:
    "envs/gatkcondaenv.yml"

4. Create a text file  *bamlist* which lists all BAM files in *input_data_dir* (see `config.yaml` and `rules/common.smk`).

5. Set the required parameters in `config.yaml`.

Note that for correct operations of such tools as VQSR (`rules/filtering.smk`) some database resources that depend on the reference genome are required.

## Train/test split (`rules/preprocessing.smk`)

All input BAMs should be labelled as belonging to either the TRAIN  or TEST set. This labelling can be provided in a csv file (*train_test_split_csv* parameter in `config.yaml`). If *train_test_split_csv* is not provided, then *N_test_bams* randomly chosen samples will be labelled as TEST samples and the remaining samples will be labelled as TRAIN samples (see `common.smk`).

As far as TEST samples are concerned, calling is performed on the entire genome. For TRAIN samples, calling is performed differently. We first generate a set of intervals that span the entire genome (see `rules/preprocessing.smk`). Then, we randomly assign each interval to one of the TRAIN samples. For each TRAIN sample, we perform calling on the corresponding intervals.

When calling is to be done on intervals generated before, use the *intervals_dir* parameter in `config.yaml`.

## Running on cluster

To run the pipeline on a SLURM cluster, run ./run_slurm.sh.
Make sure that the chosen cluster nodes have AVX support (see below).

## AVX support

Modern versions of variant callers from GATK operate much faster with AVX support.
If AVX instructions aren't supported by the CPU, calling may take ages.
Even if the CPU supports AVX instructions, make sure that the Mutect2 log file
has the message "Using CPU-supported AVX-512 instructions", otherwise
calling will be slower.
