#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_pileup
#PBS -r n

#Useful variables
ref_fasta_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta
sam_ore_bams=../bam_alignments/sam.bam\ ../bam_alignments/ore.bam
bc_bams=../bam_alignments/bc_long.bam\ ../bam_alignments/bc_short_1.bam\ ../bam_alignments/bc_short_2.bam\ ../bam_alignments/bc_short_3.bam\ ../bam_alignments/bc_short_4.bam
bsa_bams=../bam_alignments/bsa_high_1.bam\ ../bam_alignments/bsa_high_2.bam\ ../bam_alignments/bsa_low_3.bam\ ../bam_alignments/bsa_low_4.bam
vcf_file=../snp_data/all_snps.vcf
snp_list_file=snp_list.txt
pileup_file=../snp_data/all_pileup.txt

#Change the current working directory
cd $PBS_O_WORKDIR

#Generate a list of SNPs we want to output pileups for
grep '^[^#]' $vcf_file | cut -f 1,2 > $snp_list_file

#Generate the pileup
samtools mpileup -l $snp_list_file -f $ref_fasta_file $sam_ore_bams $bc_bams $bsa_bams > $pileup_file
