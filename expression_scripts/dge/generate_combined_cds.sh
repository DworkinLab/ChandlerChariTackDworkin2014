#!/bin/sh -login
#PBS -o joberrors_combined
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_combined_cds
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

sam_fasta=../resequenced_cds/sam_cds.fasta
ore_fasta=../resequenced_cds/ore_cds.fasta
output_file=../resequenced_cds/combined_cds.fasta

module load use.cus python
python generate_combined_cds.py $sam_fasta SAM $ore_fasta ORE $output_file