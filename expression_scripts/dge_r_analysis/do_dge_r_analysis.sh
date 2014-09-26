#!/bin/sh -login
#PBS -o joberrors_dge_mutant_nonmutant_analysis
#PBS -j oe
#PBS -l nodes=1:ppn=2,walltime=1:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N mutant_nonmutant_analysis
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

module load R

R --no-save < mutant_nonmutant_tests.R
R --no-save < allele_specific_tests.R
