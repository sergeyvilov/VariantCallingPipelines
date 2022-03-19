rule mutect2_pon:
    '''
    Call on normal samples only to prepare VCFs for the panel of normals (PON)
    '''
    input:
        normal_bam = input_data_dir + "{patient}.normal.bam",
        ref = refgen_fa,
        germline_resource = config["gnomAD_light_vcf"],
    output:
        vcf = pon_dir + "called/{patient}.{contig}.vcf.gz",
        tbi = pon_dir + "called/{patient}.{contig}.vcf.gz.tbi",
    params:
        max_population_af = 0.0005, #variants with lower gnomAD AF can be somatic
        normal_bam_name = lambda wildcards: bam_names[wildcards.patient + '.normal.bam'],
        java_options='-Xmx4g -Xms4g -Djava.io.tmpdir=/storage/groups/epigenereg01/workspace/projects/vale/tmp', #Setting a tmp dir solves AVX error!
    threads: 4
    log:
        logs_dir + 'mutect2_pon/{patient}/{contig}.log'
    shell:
        """
        gatk  --java-options "{params.java_options}" Mutect2 \
            -I {input.normal_bam} \
            -R {input.ref} \
            -tumor {params.normal_bam_name} \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --germline-resource {input.germline_resource} \
            --max-population-af {params.max_population_af} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} \
            -L {wildcards.contig} \
            > {log} 2>&1
        """

rule called_concat:
    '''
    Concatenate calls from different chromosomes of the same sample
    '''
    input:
        vcf = expand(pon_dir + "called/{{patient}}.{contig}.vcf.gz",contig=ref_contigs),
        tbi = expand(pon_dir + "called/{{patient}}.{contig}.vcf.gz.tbi",contig=ref_contigs),
    output:
        vcf = pon_dir + "called_concat/{patient}.vcf.gz",
        tbi = pon_dir + "called_concat/{patient}.vcf.gz.tbi",
    shell:"""
    bcftools concat {input.vcf} -Oz -o {output.vcf}; \
    tabix {output.vcf}
    """

rule gather_variants:
    '''
    Prepare a signle VCF for all the patients to create a pon
    '''
    input:
        vcfs = expand(pon_dir + "called_concat/{patient}.vcf.gz", patient=all_patients),
    output:
        vcf = pon_dir + "called_concat/merged.vcf.gz",
        tbi = pon_dir + "called_concat/merged.vcf.gz.tbi",
    shell:
        r'''
        bcftools merge {input.vcfs} -Oz -o {output.vcf};
        tabix {output.vcf}
        '''

rule create_pon:
    input:
        vcf = pon_dir + "called_concat/merged.vcf.gz",
        ref = refgen_fa
    output:
        vcf = pon_dir + "pon.vcf.gz",
        tbi = pon_dir + "pon.vcf.gz.tbi"
    log:
        logs_dir + 'CreateSomaticPanelOfNormals/log'
    shell:
        """
        gatk CreateSomaticPanelOfNormals -R {input.ref} -V {input.vcf} -O {output.vcf} > {log} 2>&1
        """
