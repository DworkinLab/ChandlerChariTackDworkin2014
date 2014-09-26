#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N make_pileup
#PBS -r n

#Useful variables
ref_fasta_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta
sam_ore_bams=../bam_alignments/sam.bam\ ../bam_alignments/ore.bam
bcf_file=../snp_data/var.raw.1.bcf
vcf_file=../snp_data/all_snps.1.vcf

#Change the current working directory
cd $PBS_O_WORKDIR

samtools mpileup -Df $ref_fasta_file -l cut.bed $sam_ore_bams > cut.pileup.txt
samtools mpileup -Df $ref_fasta_file -l sal.bed $sam_ore_bams > sal.pileup.txt
