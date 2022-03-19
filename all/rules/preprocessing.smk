rule generate_beds:
    '''
    Generate BED files defining intervals on which calling should be performed
    '''
    input:
        #TRAIN BAMS: we split the genome length into equally-sized windows and randomly assign each window to one of the train BAMs
        #TEST BAMS: calling is performed on the entire genome
        train_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = train_bams),
        test_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = test_bams),
    output:
        expand(intervals_dir + '{bam_file_name}/{chrom}.bed', bam_file_name = all_bams, chrom = chroms_GRCh37)
    params:
        train_bams = train_bams,
        test_bams = test_bams,
        workdir = intervals_dir,
        window_size = 500000, #should be small enough s.t. the composite genome contains parts of all the train samples and large enough to minimize the number of step_size-window_size gaps
        step_size = 500000 + 150, #window_size + read length+1, to avoid that pieces from different train samples overlap
        genome = 'GRCh37',
        split_by_chromosome = True,
    script:
        "../scripts/make_beds.py"

#rule index_bam:
#    input:
#        bam = input_data_dir + "{bam_file_name}.bam"
#    output:
#        bai = input_data_dir + "{bam_file_name}.bam.bai"
#    shell:
#        r'''
#        samtools index {input.bam}
#        '''
