#!/bin/sh -login
#PBS -o joberrors_combine_quantitative_data
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N combine_quantitative_data
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python combine_quantitative_data_allele_nonspecific.py ../dge_data/dge_data_all.csv ../flybase/flybase_map_table.csv ../dge_data/fc3nonspecific_counts.csv ../dge_data/fc4nonspecific_counts.csv ../dge_data/fc5nonspecific_counts.csv ../dge_data/fc6nonspecific_counts.csv
python combine_quantitative_data_allele_specific.py ../dge_data/dge_data_allele_specific_prefilter.csv ../flybase/flybase_map_table.csv ../dge_data/fc3allelespecific_counts.csv ../dge_data/fc4allelespecific_counts.csv ../dge_data/fc5allelespecific_counts.csv ../dge_data/fc6allelespecific_counts.csv