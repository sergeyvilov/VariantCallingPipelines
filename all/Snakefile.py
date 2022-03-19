include: "rules/common.smk" #define directories, make/load train/test split

rule all:
    input:
        expand(progress_dir + "called/{bam_file_name}/{chrom}.vcf.gz", bam_file_name = all_bams, chrom = chroms_GRCh37),
        expand(progress_dir + 'called_concat/{bam_file_name}.vcf.gz', bam_file_name = all_bams), #VCF files after calling
        expand(progress_dir + 'per_sample/{bam_file_name}.vcf.gz', bam_file_name = test_bams), #VCF files for each TEST sample, calling on complete genome
        progress_dir + "results/train.vcf.gz",#single VCF file for all TRAIN samples, on each sample calling is performed only on a fixed set of regions
        #progress_dir + "results/test.vcf.gz",



include: "rules/preprocessing.smk"
include: "rules/calling.smk"
include: "rules/postprocessing.smk"
