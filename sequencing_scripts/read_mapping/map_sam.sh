#!/bin/sh -login
#PBS -o joberrors_map_sam
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N map_sam
#PBS -r n

ref_fasta_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta
ref_fasta_index_file=../reseq_data/reference/dmel-all-chromosome-r5.41.fasta.fai
fastq_file_1=../reseq_data/sam_ore/SAM_1.fastq
fastq_file_2=../reseq_data/sam_ore/SAM_2.fastq
sai_file_1=../temp_alignment_data/sam_1.sai
sai_file_2=../temp_alignment_data/sam_2.sai
sam_file=../temp_alignment_data/sam.sam
unsorted_bam_file=../temp_alignment_data/sam_unsorted.bam
sorted_bam_file_no_ext=../bam_alignments/sam
sorted_bam_file=../bam_alignments/sam.bam

#Change the working directory
cd $PBS_O_WORKDIR

#Align the forward and reverse strands to the reference
bwa aln -n 4 $ref_fasta_file $fastq_file_1 > $sai_file_1
bwa aln -n 4 $ref_fasta_file $fastq_file_2 > $sai_file_2

#Generate a SAM file from the bwa alignments
bwa sampe $ref_fasta_file $sai_file_1 $sai_file_2 $fastq_file_1 $fastq_file_2 > $sam_file

#Convert the SAM file to a BAM file
samtools import $ref_fasta_index_file $sam_file $unsorted_bam_file

#Sort the BAM file
samtools sort $unsorted_bam_file $sorted_bam_file_no_ext

#Index the BAM file
samtools index $sorted_bam_file