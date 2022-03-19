#Here BQSR is implemented
#https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-
#https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR
#ATTENTION: for each sample, minimal amount of data is 100M bases (1/30 of a genome), i.e. the interval files should include enough contigs



rule recalibrate_base_qualities:
    #create recalibration tables. This doesn't create new BAM files!
    #we don't use the whole sample to create a table: this would take a lot of time
    #we rely instead on intervals randomly scattered across the genome
    input:
        bam = input_data_dir + '{sample_name}.bam',
        bai = input_data_dir + '{sample_name}.bai',
        ref = refgen_fa,
    output:
        bqsr_table = protected(recal_dir + "tables/{sample_name}.recal_table"),
    log:
        logs_dir + "gatk/bqsr/{sample_name}.log",
    params:
        java_options='--java-options  "-Xmx4g -Xms4g"',
        #we can include all resources but using only dbsnp file should work pretty well
#        known="--known-sites " + g1k_file + " --known-sites " + mills_file + \
#        " --known-sites " + knownindels_file + " --known-sites " + dbsnp_common,
        known = "--known-sites " + dbsnp_all,
        intervals=(
            '-L ' + config['params']['bqsr']['intervals_list']
            if config['params']['bqsr']['use_intervals']
            else ''
        ) #with -XL option specify intervals to exclude, with -L option specify intervals to include
    shell:
        r'''gatk {params.java_options} \
        BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        {params.known} \
        {params.intervals} \
        -O {output} > {log} 2>&1'''

#Apply BQSR to create BAM files with refined quality scores using recalibration tables
rule apply_bqsr:
    input:
        bam = input_data_dir + '{sample_name}.bam',
        table = recal_dir + "tables/{sample_name}.recal_table",
        ref = refgen_fa,
    output:
        bam = protected(recal_dir + "{sample_name}.bam"),
    log:
        logs_dir + "gatk/apply_bqsr/{sample_name}.log",
    params:
        java_options = '--java-options  "-Xmx4g -Xms4g"',
        intervals = ['-L '+contig+' ' for contig in ref_contigs], #don't apply BQSR on junk staff,e.g. decoy sequences
    shell:
        r'''gatk {params.java_options} \
        ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        {params.intervals} \
        --bqsr-recal-file {input.table} \
        -O {output} > {log} 2>&1'''

#index recalibrated BAM files
rule index_recalibrated_bam:
    input:
        bam = recal_dir + '{sample_name}.bam'
    output:
        bai = recal_dir + '{sample_name}.bai'
    log:
        logs_dir + "samtools/index_recalibrated_bam/{sample_name}.log",
    wrapper:
        "0.59.2/bio/samtools/index"
