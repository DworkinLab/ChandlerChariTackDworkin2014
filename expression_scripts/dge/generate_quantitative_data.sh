#!/bin/sh -login
#PBS -o joberrors_generate_quantitative_data
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_quantitative_data
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python generate_quantitative_data.py ../dge_data/fc3_preprocessed.csv ../dge_data/FC3CLEAN.txt ../dge_data/fc3nonspecific_counts.csv ../dge_data/fc3allelespecific_counts.csv
python generate_quantitative_data.py ../dge_data/fc4_preprocessed.csv ../dge_data/FC4CLEAN.txt ../dge_data/fc4nonspecific_counts.csv ../dge_data/fc4allelespecific_counts.csv
python generate_quantitative_data.py ../dge_data/fc5_preprocessed.csv ../dge_data/FC5CLEAN.txt ../dge_data/fc5nonspecific_counts.csv ../dge_data/fc5allelespecific_counts.csv
python generate_quantitative_data.py ../dge_data/fc6_preprocessed.csv ../dge_data/FC6CLEAN.txt ../dge_data/fc6nonspecific_counts.csv ../dge_data/fc6allelespecific_counts.csv
