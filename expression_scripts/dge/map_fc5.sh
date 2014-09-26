#!/bin/sh -login
#PBS -o joberrors_map_fc5
#PBS -j oe
#PBS -l nodes=1:ppn=2,walltime=8:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N map_fc5
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

ref_fasta_file=../resequenced_cds/combined_cds.fasta

input_fasta_file=../dge_data/fc5.fasta
sai_file=../dge_temp_alignment/fc5.sai
sam_file=../dge_data/fc5.sam

#Align to the reference CDS
bwa aln -n 2 -k 0 -l 4 -t 2 $ref_fasta_file $input_fasta_file > $sai_file
#Generate a SAM file from the bwa alignment
bwa samse -n 10000 $ref_fasta_file $sai_file $input_fasta_file > $sam_file
