#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N generate_expanded-resequenced_cds.y
#PBS -r n

cd $PBS_O_WORKDIR

module load use.cus python

python generate_expanded_resequenced_cds.py ../../expression_scripts/flybase/flybase_map_table.csv ../consensus_sequences/sam_consensus.fa ../reseq_data/reference/dmel-all-chromosome-r5.41.fasta sam_binding_seqs.fasta
python generate_expanded_resequenced_cds.py ../../expression_scripts/flybase/flybase_map_table.csv ../consensus_sequences/ore_consensus.fa ../reseq_data/reference/dmel-all-chromosome-r5.41.fasta ore_binding_seqs.fasta