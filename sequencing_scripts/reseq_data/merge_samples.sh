#!/bin/sh -login
#PBS -o joberrors
#PBS -j oe
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M cholden23@gmail.com
#PBS -m abe
#PBS -N merge_samples
#PBS -r n

#Change the working directory
cd $PBS_O_WORKDIR

gunzip -r .

cd sam_ore
cat 62G05AAXX_6_1_pf.fastq 62G2NAAXX_8_1_pf.fastq > SAM_1.fastq
cat 62G05AAXX_6_2_pf.fastq 62G2NAAXX_8_2_pf.fastq > SAM_2.fastq
cat 62G05AAXX_7_1_pf.fastq 632BKAAXX_5_1_pf.fastq > ORE_1.fastq
cat 62G05AAXX_7_2_pf.fastq 632BKAAXX_5_2_pf.fastq > ORE_2.fastq

cd ../bc

cat 62G05AAXX_8_1_BC_long_11_pf.fastq 62FYLAAXX_8_1_BC_long_11_pf.fastq > BC_long_1.fastq
cat 62G05AAXX_8_2_BC_long_11_pf.fastq 62FYLAAXX_8_2_BC_long_11_pf.fastq > BC_long_2.fastq

cat 62G05AAXX_8_1_BC_short_1_pf.fastq 62FYLAAXX_8_1_BC_short_1_pf.fastq > BC_short1_1.fastq
cat 62G05AAXX_8_2_BC_short_1_pf.fastq 62FYLAAXX_8_2_BC_short_1_pf.fastq > BC_short1_2.fastq

cat 62G05AAXX_8_1_BC_short_2_pf.fastq 62FYLAAXX_8_1_BC_short_2_pf.fastq > BC_short2_1.fastq
cat 62G05AAXX_8_2_BC_short_2_pf.fastq 62FYLAAXX_8_2_BC_short_2_pf.fastq > BC_short2_2.fastq

cat 62G05AAXX_8_1_BC_short_3_pf.fastq 62FYLAAXX_8_1_BC_short_3_pf.fastq > BC_short3_1.fastq
cat 62G05AAXX_8_2_BC_short_3_pf.fastq 62FYLAAXX_8_2_BC_short_3_pf.fastq > BC_short3_2.fastq

cat 62G05AAXX_8_1_BC_short_4_pf.fastq 62FYLAAXX_8_1_BC_short_4_pf.fastq > BC_short4_1.fastq
cat 62G05AAXX_8_2_BC_short_4_pf.fastq 62FYLAAXX_8_2_BC_short_4_pf.fastq > BC_short4_2.fastq
