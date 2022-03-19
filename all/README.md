## Call candidate variants with GATK

This is a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for preparing training, evaluation or inference data for deepSNP.

The pipeline generates VCF files with variants for given BAM samples.

Input samples should be prepared according to [GATK best practices of data-preprocessing for variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)

The pipeline uses the [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) caller in tumor-only mode, so it yields all candidate variants, including true somatic and true germline variants as well as artefacts.

Somatic variants can then be filtered out and germline variants can be labelled in output VCFs throughout [postprocessing](rules/postprocessing.smk).

See [rules](rules/) for a more detailed description of what's going on.

## Pipeline setup

1. Install  [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Create and activate the conda environment

```
conda env create -f environment.yml
conda activate vcalling
```

3. Set the required parameters in `config.yaml` (see [rules](rules/) for additional details).

## Running on cluster

To run the pipeline on a SLURM cluster, run ./run_slurm.sh.
Make sure that the chosen cluster nodes have AVX support (see below).

## AVX support

Modern versions of variant callers from GATK operate much faster with AVX support.
If AVX instructions aren't supported by the CPU, calling may take ages.
Even if the CPU supports AVX instructions, make sure that the Mutect2 log file
has the message "Using CPU-supported AVX-512 instructions", otherwise
calling will be slower.
