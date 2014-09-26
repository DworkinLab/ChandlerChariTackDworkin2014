'''
Reads one of David Tack's pre-processed FCXCLEAN.txt files and makes a fasta out of it so that
it can be used by bwa

Usage:
%% python format_sequence_reads.py <FCFile.txt> <output.fasta>

'''

import sys

input_handle = open(sys.argv[1], 'r')
output_handle = open(sys.argv[2], 'w')

first_line = True
for cur_line in input_handle:
    if not first_line:
        output_handle.write("\n")
    else:
        first_line = False
    fields = cur_line.split(", ")
    seq_tag = 'CATG' + fields[0]
    count = fields[1]
    output_handle.write(">%s\n" % seq_tag)
    output_handle.write(seq_tag)

input_handle.close()
output_handle.close()
