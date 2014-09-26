'''
Reads the allele-specific CSV files output by generate_quantitative_data.py and
combines them all into one big CSV file for statistical analysis in R

Output format:
gene, fbgn, gene_transcript_start, fc3_ore_count, fc3_sam_count, fc4_ore_count, fc4_sam_count, fc5_ore_count, fc5_sam_count, fc6_ore_count, fc6_sam_count

Usage:
%% python combine_quantitative_data_allele_specific.py <output.csv> <flybase_key.csv> <fc3nonspecific.csv> <fc4nonspecific.csv> ...

'''

import sys


###############################################
#Class to store flybase information on a gene
###############################################
class flybase_record:
    def __init__(self, fbgn="", symbol="", chrom="", start=-1, stop=-1):
        self.fbgn = fbgn
        self.symbol = symbol
        self.chrom = chrom
        self.start = start
        self.stop = stop

###############################################
#Class to store all our data
###############################################
#Next we need to do the combining...
class gene_transcript_start_record:
    def __init__(self, fbgn="", gene="", gene_transcript_start="", fc3_ore_count=0, fc4_ore_count=0, fc5_ore_count=0, fc6_ore_count=0, fc3_sam_count=0, fc4_sam_count=0, fc5_sam_count=0, fc6_sam_count=0):
        self.fbgn = fbgn
        self.gene = gene
        self.gene_transcript_start = gene_transcript_start
        self.fc3_ore_count = fc3_ore_count
        self.fc4_ore_count = fc4_ore_count
        self.fc5_ore_count = fc5_ore_count
        self.fc6_ore_count = fc6_ore_count
        self.fc3_sam_count = fc3_sam_count
        self.fc4_sam_count = fc4_sam_count
        self.fc5_sam_count = fc5_sam_count
        self.fc6_sam_count = fc6_sam_count

    def get_output_string(self):
        return ",".join([self.gene, self.fbgn, self.gene_transcript_start, str(self.fc3_ore_count), str(self.fc3_sam_count), str(self.fc4_ore_count), str(self.fc4_sam_count), str(self.fc5_ore_count), str(self.fc5_sam_count), str(self.fc6_ore_count), str(self.fc6_sam_count)])

###############################################
# Script body
###############################################

#Read in the flybase information
flybase_handle = open(sys.argv[2], 'r')
flybase = {}
for cur_line in flybase_handle:
    cur_line = cur_line.strip()
    fields = cur_line.split(",")
    cur_record = flybase_record(fields[0], fields[1], fields[2], int(fields[3]), int(fields[4]))
    flybase[fields[1]] = cur_record
flybase_handle.close()
#Keys = gene_transcript_start
#Values = flybase records

def process_fc_file(handle, which_count):
    for cur_line in handle:
        cur_line = cur_line.strip()
        fields = cur_line.split(",")
        gene = fields[0]
        if not (gene in flybase):
            print "Gene %s not found in flybase!" % gene
            continue
        gene_transcript_start = fields[1]
        ore_count = int(fields[2])
        sam_count = int(fields[3])
        fbgn = flybase[gene].fbgn
        #If the current gene is not already in the dictionary, we need to add it
        if not (gene_transcript_start in all_data):
            all_data[gene_transcript_start] = gene_transcript_start_record(fbgn, gene, gene_transcript_start, 0, 0, 0, 0, 0, 0, 0, 0)
        if which_count == 3:
            all_data[gene_transcript_start].fc3_ore_count = ore_count
            all_data[gene_transcript_start].fc3_sam_count = sam_count
        elif which_count == 4:
            all_data[gene_transcript_start].fc4_ore_count = ore_count
            all_data[gene_transcript_start].fc4_sam_count = sam_count
        elif which_count == 5:
            all_data[gene_transcript_start].fc5_ore_count = ore_count
            all_data[gene_transcript_start].fc5_sam_count = sam_count
        elif which_count == 6:
            all_data[gene_transcript_start].fc6_ore_count = ore_count
            all_data[gene_transcript_start].fc6_sam_count = sam_count
            
all_data = {} #Dictionary to hold our data
              #Keys = gene symbols
              #Values = gene_transcript_count_records

#Loop through all the flow cell files and add their information to our big dictionary of data
fc_nums = [3,4,5,6]
fc_files = sys.argv[3:7]
for (num, file) in zip(fc_nums, fc_files):
    cur_handle = open(file, 'r')
    process_fc_file(cur_handle, num)
    cur_handle.close()
    
#Now save the results
output_handle = open(sys.argv[1], 'w')
first_line = True
for gene_transcript_start, gene_info in all_data.items():
    if first_line:
        first_line = False
    else:
        output_handle.write("\n")
    output_handle.write(gene_info.get_output_string())
output_handle.close()
