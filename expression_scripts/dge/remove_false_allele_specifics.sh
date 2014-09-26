#!/bin/sh -login
#PBS -o joberrors_remove_false_allele_specifics
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=1gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N remove_false_allele_specifics
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python remove_false_allele_specifics.py ../dge_data/dge_data_allele_specific_prefilter.csv ../resequenced_cds/ore_cds.fasta ../resequenced_cds/sam_cds.fasta > ../dge_data/dge_data_allele_specific.csv