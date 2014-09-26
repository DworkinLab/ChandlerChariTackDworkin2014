#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=4gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_motifscan_pvals
#PBS -r n

cd $PBS_O_WORKDIR

module load Python
python calculate_base_comp.py sam_binding_seqs.fasta sam_gc.txt
python calculate_base_comp.py ore_binding_seqs.fasta ore_gc.txt

module load R
R --no-save < generate_motifscan_pvals.R > generate_motifscan_pvals.output.txt
