#!/bin/sh -login
#PBS -o joberrors_sam
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_resequenced_cds_sam
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

resequence_genome=../../sequencing_scripts/consensus_sequences/sam_consensus.fa
cds_file=../flybase/dmel-all-CDS-r5.41.fasta
output_file=../resequenced_cds/sam_cds.fasta

module load use.cus python
python generate_resequenced_cds.py $cds_file $resequence_genome $output_file