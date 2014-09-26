#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=20:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N index_cds
#PBS -r n

cd $PBS_O_WORKDIR

#Make a bwa index for the reference genome
bwa index ../resequenced_cds/combined_cds.fasta

#Make a samtools index of the reference genome
samtools faidx ../resequenced_cds/combined_cds.fasta
