include: "rules/common.smk"



##### Target rules #####


rule all:
    input:
        filtered_dir + "by_sample/splitted.ok"

        #vcf = progress_dir + "merged.vcf.gz",
        #tbi = progress_dir + "merged.vcf.gz.tbi",
        #expand(progress_dir + "{vartype}.vcf.gz", vartype=['snvs','indels']),

        #expand(filtered_dir + "{vartype}.recalibrated.vcf.gz", vartype=['snvs','indels']),
        #expand(progress_dir + "called_concat/{bam_file_name}.vcf", bam_file_name = train_bams),
        #expand(progress_dir + 'intervals/{bam_file_name}/{chrom}.bed', bam_file_name = all_bams, chrom = chroms_GRCh37),


include: "rules/preprocessing.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/postprocessing.smk"
