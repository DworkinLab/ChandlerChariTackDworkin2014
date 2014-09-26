#!/bin/sh -login
#PBS -o overlap_tests
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N overlap_tests
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load R

R --no-save < overlap_tests.R > overlap_tests.results.txt

