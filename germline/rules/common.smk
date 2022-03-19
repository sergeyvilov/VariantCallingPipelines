from snakemake.utils import validate
from snakemake.utils import min_version

import numpy as np
import pandas as pd
import os
import json

min_version("5.18.0")

configfile: "config.yaml"

wildcard_constraints:
    vartype="snvs|indels"

gatk_path = config["gatk_path"]

input_data_dir = config["input_data_dir"]
progress_dir = config["progress_dir"]
logs_dir = progress_dir + 'logs/'

refgen_name = config['refgen_name']
resources = config['resources'][refgen_name]
refgen_fa = resources['refgen_fa']

#create output dirs depending on the chosen filtering type
if not config["filtering"]["type"]=="vqsr":
    filtered_dir = progress_dir + 'filtered/hard/'
else:
    filtered_dir = progress_dir + f'filtered/vqsr/{config["filtering"]["vqsr"]["truth_sensitivity_filter_level"]}/'

if config["intervals_dir"]=='None':
    intervals_dir = os.path.join(progress_dir,'intervals/')
else:
    intervals_dir = config["intervals_dir"]

if refgen_name == 'GRCh37':
    chroms = [str(chr_number) for chr_number in range(1,23)]+['X', 'Y'] #GRCh37 contigs on which calling should be done

#all_bams.remove('p_0_54.tumor') #MLL project: don't call variants for this sample (according to Strelka results it gives a lot of false positives)

#we split BAM files into two groups: TRAIN and TEST BAMs
#TRAIN BAMs:
#first, we generate a set if intervals that span an entire genome
#then, we randomly assign each interval to one of the TRAIN BAMs
#for each TRAIN BAM, we perform calling on the corresponding intervals
#TEST BAMs:
#for each test BAM, we perform calling on the entire genome

if config["train_test_split_csv"]=='None':

    #if the train/test split isn't provided

    #get list of all BAMs

    with open(config['bam_list'],'r') as f:
        all_bams = f.read().split('\n')[:-1]

    #remove .bam extension
    all_bams = [bam_name.replace('.bam','') for bam_name in all_bams]

    np.random.seed(1)
    np.random.shuffle(all_bams)

    #after shuffling, allocate the first N_test_bams to the test set and the rest to the train set
    test_bams = all_bams[:config['N_test_bams']]
    train_bams = all_bams[config['N_test_bams']:]

else:

    #if the train/test split is provided

    all_bams = pd.read_csv(config["train_test_split_csv"])
    train_bams = all_bams[all_bams.split=='TRAIN']['bam'].values
    test_bams = all_bams[all_bams.split=='TEST']['bam'].values
    all_bams = all_bams['bam'].values

os.makedirs(progress_dir, exist_ok=True)

#save train/test split in the progress_dir
with open(os.path.join(progress_dir, 'train_test_split.csv'), 'w') as f:
    f.writelines("bam,split\n")
    f.writelines([bam_name+',TRAIN\n' if bam_name in train_bams else bam_name+',TEST\n' for bam_name in all_bams ])

#save configuration parameters in the progress_dir
with open(os.path.join(progress_dir, 'params.json'), 'w') as f:
    json.dump(config, f)
