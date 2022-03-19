## Train/test split (`preprocessing.smk`)


All input BAMs should be labelled as belonging to either the TRAIN (for deepSNP training) or TEST (for deepSNP evaluation or inference) set. This labelling can be provided in a csv file (*train_test_split_csv* parameter in `config.yaml`). If *train_test_split_csv* is not provided, then *N_test_bams* randomly chosen samples will be labelled as TEST samples and the remaining samples will be labelled as TRAIN samples (see `common.smk`).

As it is often important to perform evaluation/inference on all variants for a TEST sample,  we perform calling on the entire genome for each TEST sample. For TRAIN samples, calling is performed differently. Although variants from all TRAIN samples should be represented in the train set to obtain a good generalisation, calling on the whole genome for each TRAIN sample is redundant. So, we first generate a set if intervals that span the entire genome (see `preprocessing.smk`). Then, we randomly assign each interval to one of the TRAIN samples. For each TRAIN sample, we perform calling on the corresponding intervals.

When calling should be done on intervals generated before, use the *intervals_dir* parameter in `config.yaml`.

## Calling (`calling.smk`)

We call all variants using the [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) caller in tumor-only mode, so it yields all candidate variants, including true somatic and true germline variants as well as artefacts.

Parallelisation is achieved by launching Mutect2 for each chromosome separately.

## Postprocessing (`postprocessing.smk`)

Postprocessing comprises the following steps:

1. Per-chromosome VCF files are concatenated for each sample.
2. Somatic variants are removed if *somatic_vcf* is provided in  `common.smk`.
3. Only bialleic variants are selected.
4. Variants are annotated with gnomAD allele frequency (AF) if *gnomAD_light_vcf* is provided.
*gnomAD_light_vcf* should represent a lighweight version of the gnomAD database after applying gnomAD filtering,
accumulating all mutations in a single file and removing all information from the INFO field except the population allele frequency (AF).
5. Germline variants are marked if germline VCF files are provided in *germline_vcfs_dir* (note NAME PATTERN in the *mark_germline* rule).
6. Variant allele fraction (VAF) and read depth (DP) are added to the INFO field.
7. *Project* name and BAM name are added to the INFO field. Since at this stage each VCF file corresponds to one sample,
no conflicts are expected.
8. All TRAIN calls are merged to a single file. Since calling intervals of different samples don't overlap, each VCF record corresponds to a single sample, so no again no conflicts for the BAM record are expected.

Steps 3,5,6 are optional and not necessary for deepSNP data preparation.
Step 7 can also be omitted but *bam_matching_csv* should be provided when calling `generate_tensors.py`.
Step 2 can be omitted when the number of SNP variants per sample is negligible compared to the number of germline variants and artefacts (i.e. in most cases).
