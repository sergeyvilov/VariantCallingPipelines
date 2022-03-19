#!/bin/bash

#Run the variant calling pipeline on a SLURM cluster

#use -x to exclude servers with no AVX support (otherwise calling takes ages)!
#The message "Using CPU-supported AVX-512 instructions" should appear in the Mutect2 log
#If this message is absent, the node should also be excluded as AVX works somehow more slowly on these machines

source /home/icb/sergey.vilov/.bashrc
conda activate vale-bio
rm -f slurm_logs/*
mkdir slurm_logs
snakemake -s Snakefile.py -k --use-conda --latency-wait 180 --cluster-config cluster.yaml --cluster 'sbatch   -p cpu_p \
-x ibis-ceph-0[17-19],ibis-ceph-0[02-05,08-16],ibis216-010-0[35,37,64],ibis216-010-0[51],ibis216-010-0[68-70],ibis216-010-0[71],ibis216-224-0[10-11] \
--mem={cluster.mem} --threads-per-core={cluster.threads_per_core} --nice=10000 -c {cluster.cores} -o slurm_logs/%j.out' -j 300
