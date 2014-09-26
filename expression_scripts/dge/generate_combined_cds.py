'''
Merges two fasta files into one, adding a prefix to the name of each sequence
(different prefix for each of the two files)

Usage:
%% python generate_combined_cds.py <fasta1.fasta> <prefix1> <fasta2.fasta> <prefix2> <output.fasta>

Requirements: Biopython must be installed
'''

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

combined_seqs = []

def read_seqs(prefix, filename):
    seq_handle = open(filename, "rU")
    seqs = []
    for record in SeqIO.parse(seq_handle, "fasta"):
        #record.description = record.description.partition(record.id + " ")[2]
        record.description = '' #For some reason we empty out the description, or else bwa has problems later on
        record.id = "%s-%s" % (prefix, record.id)
        seqs.append(record)
    seq_handle.close()
    return seqs

#Read the first fasta file
prefix1 = sys.argv[2]
seqs1 = read_seqs(prefix1, sys.argv[1])

#Read the second fasta file
prefix2 = sys.argv[4]
seqs2 = read_seqs(prefix2, sys.argv[3])

output_handle = open(sys.argv[5], 'w')
SeqIO.write(seqs1, output_handle, "fasta")
SeqIO.write(seqs2, output_handle, "fasta")
output_handle.close()