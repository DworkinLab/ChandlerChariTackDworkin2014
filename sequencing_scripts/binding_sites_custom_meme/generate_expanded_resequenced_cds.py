'''
Uses the flybase map table, the original reference genome, and the re-sequenced genome to generate
a fasta file with gene sequences along with 10kb of upstream and downstream flanking sequences for
each gene

Usage:
%% python generate_expanded_resequenced_cds.py <flybase_map_table.fasta> <resequenced_genome.fa> <reference_genome.fasta> <output.fasta>

'''

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

#Ignore any genes that aren't in these chromosomes:
valid_chroms = ["X", "2L", "2R", "3L", "3R"]

print "Reading the reference genome..."

#Read the original reference file
ref_handle = open(sys.argv[3], "rU")
ref_seqs = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
ref_handle.close()

print "Reading the re-sequenced genome..."

#Read the re-sequenced genome file
rsq_seqs = []
rsq_handle = open(sys.argv[2], "rU")
rsq_seqs = SeqIO.to_dict(SeqIO.parse(rsq_handle, "fasta"))
rsq_handle.close()

#Replaces all the N's in a sequence with the appropriate base from the template
#Assumes it's already all uppercase
def replace_ns(new_seq, template):
    if len(new_seq) != len(template):
        if len(new_seq) < len(template):
            new_seq = new_seq[0:len(template)]
        else:
            raise Exception("replace_ns: original and template are not the same length")
    result = new_seq.upper() #str(new_seq)
    if "N" in result:
        for i in range(len(result)):
            if result[i] == "N":
                result = result[0:i] + template[i:i+1] + result[i+1:len(result)]
    return result


output_handle = open(sys.argv[4], 'w')

q = 0

print "Processing..."

#Go through the flybase map table
flybase_handle = open(sys.argv[1], 'r')
for cur_line in flybase_handle:
    q += 1
    if q % 100 == 0:
        print q
    cur_fields = cur_line.strip().split(",")
    fbgn = cur_fields[0]
    symbol = cur_fields[1]
    chrom = cur_fields[2]
    start = int(cur_fields[3]) - 10000
    if start < 0:
        start = 0
    stop = int(cur_fields[4]) + 10000
    #Pull out the sequence
    if chrom in valid_chroms:
        ref_version = ref_seqs[chrom].seq[start:stop]
        rsq_version = rsq_seqs[chrom].seq[start:stop]
        output_seq = str(replace_ns(rsq_version, ref_version))
        output_handle.write(">%s\n" % fbgn)
        output_handle.write("%s\n" % output_seq)

output_handle.close()        
flybase_handle.close()