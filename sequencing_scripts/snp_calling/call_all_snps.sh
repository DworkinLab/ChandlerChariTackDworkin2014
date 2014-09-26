#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N call_snps
#PBS -r n

#Useful variables
ref_fasta_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta
sam_ore_bams=../bam_alignments/sam.bam\ ../bam_alignments/ore.bam
bc_bams=../bam_alignments/bc_long.bam\ ../bam_alignments/bc_short_1.bam\ ../bam_alignments/bc_short_2.bam\ ../bam_alignments/bc_short_3.bam\ ../bam_alignments/bc_short_4.bam
bsa_bams=../bam_alignments/bsa_high_1.bam\ ../bam_alignments/bsa_high_2.bam\ ../bam_alignments/bsa_low_3.bam\ ../bam_alignments/bsa_low_4.bam
bcf_file=../snp_data/var.raw.bcf
vcf_file=../snp_data/all_snps.vcf

#Change the current working directory
cd $PBS_O_WORKDIR

samtools mpileup -uDf $ref_fasta_file $sam_ore_bams $bc_bams $bsa_bams | bcftools view -bvcg - > $bcf_file
bcftools view $bcf_file | vcfutils.pl varFilter -D100 > $vcf_file