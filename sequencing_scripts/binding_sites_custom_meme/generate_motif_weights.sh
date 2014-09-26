#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_motif_weights
#PBS -r n

cd $PBS_O_WORKDIR

MotifWeight sam_binding_seqs.fasta sd.wtmx -m 3 -mt 0.99 > motifweights_sam_custom_meme.txt
MotifWeight ore_binding_seqs.fasta sd.wtmx -m 3 -mt 0.99 > motifweights_ore_custom_meme.txt