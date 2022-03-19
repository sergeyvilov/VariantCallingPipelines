rule mutect2:
    '''
    Call somatic variants on a given chromosome for a given tumor-normal pair
    '''
    input:
        tumor_bam = input_data_dir + "{patient}.tumor.bam",
        normal_bam = input_data_dir + "{patient}.normal.bam",
        ref = refgen_fa,
        germline_resource = config["gnomAD_light_vcf"],
        pon = pon_dir + "pon.vcf.gz",
    output:
        vcf = called_dir + "{patient}.{contig}.vcf.gz",
        stats = called_dir + "{patient}.{contig}.vcf.gz.stats",
    params:
        normal_bam_name = lambda wildcards: bam_names[wildcards.patient + '.normal.bam'],
        af_of_alleles_not_in_resource = 6.6e-6, #1/(2*number_of_samples_in_gnomAD)
        java_options='-Xmx4g -Xms4g -Djava.io.tmpdir=/storage/groups/epigenereg01/workspace/projects/vale/tmp', #Setting a tmp dir solves AVX error!
    threads: 4
    log:
        logs_dir + 'mutect2/{patient}/{contig}.log'
    shell:
        """
        gatk  --java-options "{params.java_options}" Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -L {wildcards.contig} \
            -normal {params.normal_bam_name} \
            --panel-of-normals {input.pon} \
            --germline-resource {input.germline_resource} \
            --native-pair-hmm-threads {threads} \
            --af-of-alleles-not-in-resource {params.af_of_alleles_not_in_resource} \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            -O {output.vcf} > {log} 2>&1
        """

rule concat_calls:
    '''
    Concatenate calls from different chromosomes of the same sample
    '''
    input:
        vcf = expand(called_dir + "{{patient}}.{contig}.vcf.gz",contig=ref_contigs),
    output:
        vcf = called_dir + "{patient}.vcf.gz",
        tbi = called_dir + "{patient}.vcf.gz.tbi",
    shell:
        r'''
        bcftools concat  --no-version {input.vcf} |bcftools sort - |bgzip -c > {output.vcf};
        tabix -f {output.vcf};
        '''

rule merge_stats:
    '''
    Merge statistics from different chromosomes of the same sample
    '''
    input:
        stats = expand(called_dir + "{{patient}}.{contig}.vcf.gz.stats",contig=ref_contigs),
    output:
        stats = called_dir + "{patient}.vcf.gz.stats",
    params:
        input_stats = lambda wildcards,input: ['-stats '+stat for stat in input.stats]
    shell:
        r'''
        gatk MergeMutectStats {params.input_stats} -O {output.stats}
        '''

rule pileup_summaries:
    input:
        tumor_bam = input_data_dir + "{patient}.{sample_type}.bam",
    params:
        germline_resource = config["gnomAD_short"],
        java_options='-Xmx12g -Xms12g', #minimum 24GB for 12mln variants in gnomAD
    output:
        table = progress_dir + "pileup_summaries/{patient}.{sample_type}.pileupsummaries.table"
    log:
        logs_dir + 'pileup_summaries/{patient}.{sample_type}.log'
    shell:
        """
        gatk  --java-options "{params.java_options}" GetPileupSummaries \
            -I {input.tumor_bam} \
            -V {params.germline_resource} \
            -L {params.germline_resource} \
            -O {output.table} > {log} 2>&1
        """


rule calculate_contamination:
    input:
        tumor_table = progress_dir + "pileup_summaries/{patient}.tumor.pileupsummaries.table",
        normal_table = progress_dir + "pileup_summaries/{patient}.normal.pileupsummaries.table",
    output:
        table = progress_dir + "calculate_contamination/{patient}.contamination.table"
    log:
        logs_dir + 'calculate_contamination/{patient}.log'
    shell:
        """
        gatk CalculateContamination \
            -I {input.tumor_table}  \
            -matched {input.normal_table} \
            -O {output.table} > {log} 2>&1
        """


rule filter_calls:
    input:
        vcf = called_dir + "{patient}.vcf.gz",
        stats = called_dir + "{patient}.vcf.gz.stats",
        ref = refgen_fa,
        contamination = progress_dir + "calculate_contamination/{patient}.contamination.table"
    output:
        vcf = filtered_dir + "{patient}.filtered.vcf.gz",
        tbi = filtered_dir + "{patient}.filtered.vcf.gz.tbi",
    shell:
        """
        gatk FilterMutectCalls \
            -V {input.vcf} \
            -R {input.ref} \
            --contamination-table {input.contamination} \
            --stats {input.stats} \
            -O {output.vcf}
        """

# rule orientation_bias:
#     input:
#         get_orientationbias_input
#     output:
#         "vcfs/{patient}.read_orientation_model.tar.gz"
#     params:
#         i=lambda wildcards, input: ['-I ' + d for d in input]
#     shell:
#         """
#         gatk LearnReadOrientationModel {params.i} -O {output}
#         """
#
