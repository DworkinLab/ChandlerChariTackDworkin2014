'''
Reads in a fasta file
Outputs a tab-delimited file with the sequence name, the GC content, and the number of non-N bases

Usage:
python calculate_base_comp.py <input_file> <output_file>
'''

import sys

input_handle = open(sys.argv[1], 'r')
output_handle = open(sys.argv[2], 'w')

done = False
while not done:
    first_line = input_handle.readline().strip()
    second_line = input_handle.readline().strip()
    if (len(first_line) == 0) | (len(second_line) == 0):
        done = True
    else:
        seq_name = first_line[1:len(first_line)]
        seq = second_line.upper()
        gc_sum = seq.count('G') + seq.count('C')
        at_sum = seq.count('A') + seq.count('T')
        if (gc_sum + at_sum) < 1:
            gc_content = 0.0
        else:
            gc_content = float(gc_sum) / float(gc_sum + at_sum) #This approach makes it ignore N's, which is what we want
        output_line = '\t'.join([seq_name, str(gc_content), str(gc_sum+at_sum)])
        #print(output_line)
        output_handle.write(output_line + '\n')

input_handle.close()
output_handle.close()