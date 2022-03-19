## Variant calling pipelines

- [germline](germline/) - calling germline variants with GATK
- [somatic](somatic/) - calling somatic variants with GATK
- [all](all/) - calling somatic and germline variants as well as artefacts with GATK

These pipelines were initially designed to prepare data for deepSNP but they can also be fully or partially used
separately to achieve other goals.

All pipelines are implemented with [Snakemake](https://snakemake.readthedocs.io/en/stable/).
Each pipeline is provided with a bash script to run on a SLURM computational cluster.

Although the pipelines can also be run on a standalone machine, they make use of parallelisation, so on a computational cluster they will be more efficient.
