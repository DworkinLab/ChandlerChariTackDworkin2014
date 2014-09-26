'''
Extracts a particular chromosome/sequence from a fasta file
Usage:
python extract_chrom.py <fasta-file> <chrom>
'''

import sys

fasta_file = sys.argv[1]
chrom = sys.argv[2]

fasta_handle = open(fasta_file, 'r')

#Input the first line of the file
cur_line = fasta_handle.readline()
while (cur_line.strip() != ('>' + chrom)): #Loop through until we find the chromosome/seq we want
    cur_line = fasta_handle.readline()
    if (cur_line == ''):
        break

print(cur_line.strip()) #Print the first line of the chrom/seq
cur_line = fasta_handle.readline() #Read the next line
if cur_line != '':
    while (cur_line[0] != '>'): #Check if we've encountered a new sequence
        print(cur_line.strip()) #If we haven't, print out the current line and read the next one
        cur_line = fasta_handle.readline()
        if (cur_line == ''):
            break

fasta_handle.close() #We're all done, so close