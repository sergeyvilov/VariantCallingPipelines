## Calling somatic variants using GATK

This is a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for calling somatic variants according to [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-).

The pipeline generates VCF files with variant calls for given tumor-normal pairs.

The input consists of tumor-normal pairs, whose file names should obey the following patterns:
patient_N.normal.bam
patient_N.tumor.bam

Input samples should be prepared according to [GATK best practices of data-preprocessing for variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)

For details of the calling procedure see
https://gatk.broadinstitute.org/hc/en-us/articles/360035889791--How-to-Call-somatic-mutations-using-GATK4-Mutect2-Deprecated-
(step *FilterByOrientationBias* is omitted)

In addition, the pipeline implements quality control on BAM files is optional and can also be omitted.

After calling, the following postprocessing steps are impelmented (`rules/postprocessing.smk`):

1. Annotate calls with gnomAD population allele frequency (AF)
2. Add tumour variant allele fraction (VAF) to the INFO field
3. Remove FORMAT field
4. Add BAM name to the INFO field
5. Concatenate calls on all samples
6. Split calls into SNPs and INDELs

## Pipeline setup

1. Install  [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Create and activate the conda environment

```
conda env create -f environment.yml
conda activate vcalling
```

3. Create a text file  *all_patients* listing all *patient_N* in *input_data_dir* (see `rules/common.smk`).

4. Set the required parameters in `config.yaml`.

## Running on cluster

To run the pipeline on a SLURM cluster, run ./run_slurm.sh.
Make sure that the chosen cluster nodes have AVX support (see below).

## AVX support

Modern versions of variant callers from GATK operate much faster with AVX support.
If AVX instructions aren't supported by the CPU, calling may take ages.
Even if the CPU supports AVX instructions, make sure that the Mutect2 log file
has the message "Using CPU-supported AVX-512 instructions", otherwise
calling will be slower.
