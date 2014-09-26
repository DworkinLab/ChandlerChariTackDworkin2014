#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N motif_scan_ore_double
#PBS -r n

fasta_file=ore_binding_seqs.fasta
motif_file=sd.wtmx
output_dir=output_ore_custom_meme
motif_weight=0.0000352758

cd $PBS_O_WORKDIR

mkdir $output_dir
MotifScan $fasta_file $motif_file $motif_weight -m 3 -mt 0.99 -od $output_dir
