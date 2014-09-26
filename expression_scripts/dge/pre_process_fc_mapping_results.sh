#!/bin/sh -login
#PBS -o joberrors_preprocess_fc_mapping_results
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N preprocess_fc_mapping_results
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python pre_process_fc_mapping_results.py ../dge_data/fc3.sam ../dge_data/fc3_preprocessed.csv ../dge_data/fc3_problem_hits.csv
python pre_process_fc_mapping_results.py ../dge_data/fc4.sam ../dge_data/fc4_preprocessed.csv ../dge_data/fc4_problem_hits.csv
python pre_process_fc_mapping_results.py ../dge_data/fc5.sam ../dge_data/fc5_preprocessed.csv ../dge_data/fc5_problem_hits.csv
python pre_process_fc_mapping_results.py ../dge_data/fc6.sam ../dge_data/fc6_preprocessed.csv ../dge_data/fc6_problem_hits.csv