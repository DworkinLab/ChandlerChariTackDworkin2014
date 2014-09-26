#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N sam_consensus
#PBS -r n

#Useful variables
ref_fasta_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta
ref_fasta_index_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta.fai
input_bam_file=../bam_alignments/sam.bam
output_file=../consensus_sequences/sam_consensus.fq
output_fasta_file=../consensus_sequences/sam_consensus.fa

#Change the working directory
cd $PBS_O_WORKDIR

#Call the consensus sequence
samtools mpileup -uf $ref_fasta_index_file $input_bam_file | bcftools view -cg - | vcfutils.pl vcf2fq > $output_file

#Convert to fasta from fastq
module load use.cus python
python consensus_fq2fa.py $output_file $output_fasta_file