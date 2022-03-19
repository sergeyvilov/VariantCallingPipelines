rule mutect2:
    '''
    Call candidate variants (somatic+germline+artefacts) on selected intervals of a given tumor_bam
    '''
    input:
        tumor_bam = input_data_dir + "{bam_file_name}.bam",
        intervals = intervals_dir + "{bam_file_name}/{chrom}.bed",
        ref = refgen_fa,
        #bai = input_data_dir + "{bam_file_name}.bam.bai",
    output:
        vcf = progress_dir + "called/{bam_file_name}/{chrom}.vcf.gz",
        #stats = progress_dir + "called/{bam_file_name}/{chrom}.vcf.gz.stats",
    params:
        java_options='-Xmx4g -Xms4g -Djava.io.tmpdir=/storage/groups/epigenereg01/workspace/projects/vale/tmp', #Setting a tmp dir solves AVX error!
    log:
        progress_dir + 'logs/mutect2/{bam_file_name}/{chrom}.log'
    threads: 4 # should be equal to number of CPUs per computational node
    shell:
        """
        gatk  --java-options "{params.java_options}" Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -L {input.intervals} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} > {log} 2>&1
        """
