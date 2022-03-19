rule generate_beds:
    '''
    Generate BED files defining intervals on which calling should be performed
    '''
    input:
        #TRAIN BAMS: we split the genome length into equally-sized windows and randomly assign each window to one of the train BAMs
        #TEST BAMS: for each test BAM each chromosome will be called in its entirety i.e. intervals are the same for all test BAMs
        train_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = train_bams),
        test_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = test_bams)
    output:
        expand(intervals_dir + '{bam_file_name}/{chrom}.bed', bam_file_name = all_bams, chrom = chroms)
    params:
        train_bams = train_bams,
        test_bams = test_bams,
        workdir = intervals_dir,
        window_size = 500000, #should be small enough s.t. the composite genome is representative for all the samples and large enough to minimize the number of step_size-window_size fragments
        step_size = 500000 + 150, #window_size + read length+1, to avoid that different samples overlap
        genome = 'GRCh37',
        split_by_chromosome = True,
    script:
        "../scripts/make_beds.py"

rule bam_sample_name_table:
    '''
    Generate the correspondence table between BAM file name and sample name inside this BAM
    '''
    input:
        train_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = train_bams),
        test_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = test_bams)
    output:
        bam_names_csv = progress_dir + 'bam_names.csv',
    shell:
        r'''
        echo >{output} #creat/clear output file
        for bam in {input.train_bams}; do
            bam_sample_name=$(samtools view -H $bam|grep -oE  -m1 "SM:[^[:space:]]*"|sed 's/^SM://')
            echo  "$bam,$bam_sample_name,TRAIN" >> {output};
        done
        for bam in {input.test_bams}; do
            bam_sample_name=$(samtools view -H $bam|grep -oE  -m1 "SM:[^[:space:]]*"|sed 's/^SM://')
            echo  "$bam,$bam_sample_name,TEST" >> {output};
        done
        '''
