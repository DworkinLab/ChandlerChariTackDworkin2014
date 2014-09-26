#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=2gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N extract_all_chroms
#PBS -r n

cd $PBS_O_WORKDIR

module load use.cus python

python extract_chrom.py ../consensus_sequences/ore_consensus.fa X > ../consensus_sequences/ore_X.fa
python extract_chrom.py ../consensus_sequences/ore_consensus.fa 2L > ../consensus_sequences/ore_2L.fa
python extract_chrom.py ../consensus_sequences/ore_consensus.fa 2R > ../consensus_sequences/ore_2R.fa
python extract_chrom.py ../consensus_sequences/ore_consensus.fa 3L > ../consensus_sequences/ore_3L.fa
python extract_chrom.py ../consensus_sequences/ore_consensus.fa 3R > ../consensus_sequences/ore_3R.fa

python extract_chrom.py ../consensus_sequences/sam_consensus.fa X > ../consensus_sequences/sam_X.fa
python extract_chrom.py ../consensus_sequences/sam_consensus.fa 2L > ../consensus_sequences/sam_2L.fa
python extract_chrom.py ../consensus_sequences/sam_consensus.fa 2R > ../consensus_sequences/sam_2R.fa
python extract_chrom.py ../consensus_sequences/sam_consensus.fa 3L > ../consensus_sequences/sam_3L.fa
python extract_chrom.py ../consensus_sequences/sam_consensus.fa 3R > ../consensus_sequences/sam_3R.fa
