#Snakemake pipeline to call somatic variants

include: "rules/common.smk"

rule all:
    input:
         #qc_report = qc_dir + "multiqc_report.html",  #quality control on BAM files
         #qc_depths_table = qc_dir + "depths.csv", #quality control on BAM files
         pon = pon_dir + "pon.vcf.gz", #panel of normals
         pileupsummaries =  expand(progress_dir + "pileup_summaries/{patient}.{sample_type}.pileupsummaries.table",patient=all_patients,sample_type=['tumor','normal']),
         called =  expand(called_dir + "{patient}.vcf.gz",patient=all_patients), #Mutect calls
         filtered = expand(filtered_dir + "{patient}.filtered.vcf.gz", patient=all_patients), #filtered Mutect calls
         vcfs =  expand(progress_dir + 'results/{vartype}/{vartype}.vcf.gz', vartype=['snvs','indels']), #Concatenated vcfs from all samples, after filtering



#include: "rules/qc.smk" #quality control for BAM files
include: "rules/ref.smk" #index original BAM files
include: "rules/pon.smk" #create panel of normals
include: "rules/calling.smk" #calling, calculating contamination, filtering
include: "rules/postprocessing.smk" #postprocessing
