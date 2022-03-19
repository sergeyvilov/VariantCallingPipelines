from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.yaml"

input_data_dir = config["input_data_dir"]
progress_dir = config["progress_dir"]
recal_dir = config["input_data_dir"] + "recal/"
qc_dir = progress_dir + "qc/"
pon_dir = progress_dir + "pon/"
logs_dir = progress_dir + "logs/"
called_dir = progress_dir + "called/"
filtered_dir = progress_dir + "filtered/"

refgen_fa = config["refgen_fa"]
refgen_dict = config["refgen_dict"]

ref_contigs = [str(chr_number) for chr_number in range(1,23)]+['X', 'Y', 'MT'] #GRCh37 contigs on which calling is done

with open(input_data_dir + 'all_patients') as f:
    all_patients = f.read().split('\n')[:-1]

all_patients.remove('p_0_54') #don't call variants for this sample (according to Strelka results it gives a lot of false positives)

#get the list of BAM sample names
#to generate the list, descend to the folder with BAM files and execute the following command:
#for f in $(ls *.bam); do echo $f;samtools view -H $f|grep '@RG'|head -n 1|sed -E 's/.*SM:([^\t]*).*/\1/'; done|awk 'BEGIN {OFS="\t"} {if (NR%2==1) {sample_name=$1} else {print sample_name, $1}}'> BAM_names
with open(input_data_dir + 'BAM_names') as f:
    bam_names = f.read().split('\n')

bam_names = {bam_name.split()[0]:bam_name.split()[1] for bam_name in bam_names[:-1]} #file_name:BAM_name

sample_types = ('tumor', 'normal')
