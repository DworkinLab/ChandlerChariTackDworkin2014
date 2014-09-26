#!/bin/sh -login
#PBS -o joberrors_microarray_analysis
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N microarray_analysis
#PBS -r n

#Change the working directory
cd ${PBS_O_WORKDIR}

module load R

R --no-save < generate_illumina_table.R
R --no-save < generate_dgrc_table.R
