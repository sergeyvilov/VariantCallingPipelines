rule index_original_bam:
    '''
    Tabix original BAM files
    '''
    input:
        bam = input_data_dir + "{patient}.{sample_type}.bam",
    output:
        bai = input_data_dir + "{patient}.{sample_type}.bai",
    log:
        logs_dir + "samtools/index_original_bam/{patient}.{sample_type}.log",
    wrapper:
        "0.59.2/bio/samtools/index"
