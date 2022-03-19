def makewindows(genome,
                window_size,
                step_size=None,
                ):
    '''
    Mimics bedtools makewindows function
    window_size
		Divide each input interval (either a chromosome or a BED interval)
		to fixed-sized windows (i.e. same number of nucleotide in each window).
		Can be combined with -s <step_size>
    step_size
		Step size: i.e., how many base pairs to step before
		creating a new window. Used to create "sliding" windows.
    genome
        changes chromosome names
        genome={'GRCh37','GRCh38'}
    '''
    if step_size == None:
        step_size = window_size

    chromosome_length = {
    '1':249250621, '2':243199373, '3':198022430, '4':191154276, '5':180915260,
    '6':171115067, '7':159138663, '8':146364022, '9':141213431, '10':135534747,
    '11':135006516, '12':133851895, '13':115169878, '14':107349540, '15':102531392,
    '16':90354753, '17':81195210, '18':78077248, '19':59128983, '20':63025520,
    '21':48129895, '22':51304566, 'X':155270560, 'Y':59373566, 'MT':16569,
    } #GRCh37 contig lengths: https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37

    if genome == 'GRCh38':
        chromosome_length = {'chr'+k.replace('MT','M'):v for k,v in chromosome_length.items()} #most common GRCh38 notation

    all_intervals = []
    for chrom,chrom_length in chromosome_length.items():
        for start_pos in list(range(0,chrom_length,step_size)):
            end_pos = min(start_pos + window_size, chrom_length-1)
            if start_pos+step_size>chrom_length-1:
                end_pos=chrom_length-1
            #print(chrom,start_pos,end_pos)
            all_intervals.append({'chrom':chrom, 'start':start_pos, 'stop':end_pos})
            if end_pos==chrom_length-1:
                break

    return pd.DataFrame(all_intervals)

import pandas as pd
import numpy as np
import re
import os

def distribute_bams(bam_names, genome, window_size, step_size=None):

    '''
        Split the genome into equally-sized windowsand assign each window randomly to one of the given BAMs
    '''

    intervals = makewindows(genome, window_size, step_size) #create windows over the genome length

    N_intervals = len(intervals)
    N_bams = len(bam_names)

    #print(N_intervals)

    #each interval is to be associated with one of the BAM names
    #we don't bootstrap as we want that all sample
    #are represented at equal proportions in the composite genome

    bam_column = np.tile(bam_names, N_intervals//N_bams+1)
    bam_column = bam_column[:N_intervals]

    np.random.seed(1) #fix seed fo reproducibility
    np.random.shuffle(bam_column) #randomly assign sample to an interval

    intervals['bam_name'] = bam_column

    return intervals

def write_beds(output_dir, intervals, split_by_chromosome = False):
    '''
        For each BAM file create a BED file with the assigned intervals
    '''
    for bam_name in intervals.bam_name.unique():

        if split_by_chromosome:

            for chrom_name in intervals.chrom.unique():

                bed_dir = os.path.join(output_dir, bam_name)
                os.makedirs(bed_dir,exist_ok=True)
                bed_name = os.path.join(bed_dir, chrom_name + '.bed')

                intervals[(intervals.bam_name == bam_name) & (intervals.chrom == chrom_name)][['chrom','start','stop']].to_csv(bed_name,header=None,index=None,sep="\t")

        else:
            bed_dir = os.path.join(output_dir, re.sub('/[^/]+$','',bam_name)) #if bam name includes a directory which must be created
            os.makedirs(bed_dir,exist_ok=True)
            bed_name = os.path.join(output_dir, bam_name + '.bed')
            intervals[intervals.bam_name == bam_name][['chrom','start','stop']].to_csv(bed_name,header=None,index=None,sep="\t")

print('Running OK')

train_bams = snakemake.params.train_bams #list of BAM files which make up the train composite genome
test_bams = snakemake.params.test_bams #list of BAM files for which we will call over the entire genome
genome = snakemake.params.genome #'GRCh37' or 'GRCh38'
window_size = snakemake.params.window_size #size of each window of the composite genome
step_size = snakemake.params.step_size #how many base pairs to step before creating a new window
output_dir = snakemake.params.workdir
split_by_chromosome = snakemake.params.split_by_chromosome


# test_bams = []
# train_bams = ['A','B','C']
# genome = 'GRCh37'
# output_dir = '.'
# window_size = 100000
# step_size = 100000
# split_by_chromosome = True

# print('Test bams:',test_bams)
# print('Train bams:',train_bams)
# print('genome:',genome)
# print('output_dir:',output_dir)
# print('window_size:',window_size)
# print('step_size:',step_size)
# print('split_by_chromosome:',split_by_chromosome)

#remove path and extention
train_bams = [bam_name.replace('.bam','') for bam_name in train_bams]
test_bams = [bam_name.replace('.bam','') for bam_name in test_bams]

if train_bams and len(train_bams)>0:

    #we split the genome length into equally-sized windows and randomly assign each window to one of the train BAMs

    train_intervals = distribute_bams(train_bams, genome, window_size, step_size)

    write_beds(output_dir, train_intervals, split_by_chromosome)

if test_bams and len(test_bams)>0:

    #for each test BAM each chromosome will be called in its entirety
    #i.e. intervals are the same for all test BAMs

    test_intervals = makewindows(genome, 1000000000, 1000000000) #we basically take entire chromosomes

    for bam_name in test_bams:

        test_intervals['bam_name'] = bam_name

        write_beds(output_dir, test_intervals, split_by_chromosome)
