#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_allele_freq_data
#PBS -r n

#Useful variables
vcf_file=../snp_data/all_snps.vcf
pileup_file=../snp_data/all_pileup.txt
csv_file=../snp_data/all_allele_freq_plot_data.csv

#Change the current working directory
cd $PBS_O_WORKDIR

#Load the python module
module load use.cus python

#Generate the file
python generate_allele_freq_plot_data.py $vcf_file $pileup_file > $csv_file
