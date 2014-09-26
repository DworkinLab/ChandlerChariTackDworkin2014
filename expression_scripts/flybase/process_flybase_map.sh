#!/bin/sh -login
#PBS -o joberrors_process_flybase_map
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N process_flybase_map
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load use.cus python

python process_flybase_map.py gene_map_table_fb_2011_09.tsv flybase_map_table.csv