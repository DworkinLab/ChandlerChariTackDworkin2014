#!/bin/sh -login
#PBS -o joberrors_format_sequence_reads
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=1gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N format_sequence_reads
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python format_sequence_reads.py ../dge_data/FC3CLEAN.txt ../dge_data/fc3.fasta
python format_sequence_reads.py ../dge_data/FC4CLEAN.txt ../dge_data/fc4.fasta
python format_sequence_reads.py ../dge_data/FC5CLEAN.txt ../dge_data/fc5.fasta
python format_sequence_reads.py ../dge_data/FC6CLEAN.txt ../dge_data/fc6.fasta