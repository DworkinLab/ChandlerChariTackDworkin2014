'''
Reads a CSV file output by pre_process_fc_mapping_results.py, along with the corresponding
flow cell count file
[Which has, for a single flow cell:
  sequence_tag, gene_name, gene_transcript, start_pos, background(SAM/ORE/NA)
 ] 
And outputs two CSV files, one containing all informative transcripts (no allele specificity),
and one for allele-specific data.

The non-allele specific dataset has:
  gene_name, gene_transcript_start, count
  
The allele-specific dataset has:
  gene_name, gene_transcript_start, ORE_count, SAM_count

Usage:
%% python generate_quantitative_data.py <input_preprocessed_file.csv> <input_fc_count_file.csv> <output_nonspecific.csv> <output_allelespecific.csv>

'''

import sys

###############################################
#Class to store information about a mapping hit
#Note: Slightly different fro the one we used before!
###############################################
class mapping_hit:
    def __init__(self, tag="", gene="", transcript="", pos=-1, background=""):
        self.tag = tag
        self.gene = gene
        self.transcript = transcript
        self.pos = pos
        self.background = background
    def get_id(self):
        result = "%s-%s-%d" % (self.gene, self.transcript, self.pos)
        return result


###############################################
#The script
###############################################

#Read the FC file and create a dictionary
#   keys = sequence tags
#   values = counts
fc_count_handle = open(sys.argv[2], 'r')
fc_count_dict = {}
for cur_line in fc_count_handle:
    cur_line = cur_line.strip()
    tag = 'CATG' + cur_line.split(", ")[0]
    count = int(cur_line.split(", ")[1])
    fc_count_dict[tag] = count
fc_count_handle.close()

#Make an empty dictionary to store all the output results 
#   keys = gene-transcript-start_pos
#   values = lists of hits
#Go through the pre-processed CSV file line-by-line
#For each line
#   Add current hit to the dictionary
hit_handle = open(sys.argv[1], 'r')
hit_dict = {}
for cur_line in hit_handle:
    cur_line = cur_line.strip()
    cur_fields = cur_line.split(",")
    cur_hit = mapping_hit(cur_fields[0], cur_fields[1], cur_fields[2], int(cur_fields[3]), cur_fields[4])
    cur_id = cur_hit.get_id()
    if not (cur_id in hit_dict):
        hit_dict[cur_id] = []        
    cur_hit_list = hit_dict[cur_id]
    cur_hit_list.append(cur_hit)     
hit_handle.close()

# for cur_hit_list_key in hit_dict:
#     cur_hit_list = hit_dict[cur_hit_list_key]
#     cur_num = len(cur_hit_list)
#     if cur_num > 0:
#         print "%s has %d hits" % (cur_hit_list[0].get_id(), cur_num)
#     else:
#         print "A 0-length list..."



#Now:
#Iterate over the dictionary of hits
#For each gene-transcript-start:
#   Reset flag nas_found = False
#   Go through all the tags that mapped to it
#     Record the ORE sum
#     Record the SAM sum
#     Record the NA sum (if it's NA, also set nas_found = True)
#   Is it allele-specific? (nas_found == False?)
#   If so:
#     Write the ORE and SAM sums in the allele-specific output file
#   Write the total sum (NA + ORE + SAM) in the non-specific output file

nonspecific_handle = open(sys.argv[3], 'w')
allelespecific_handle = open(sys.argv[4], 'w')

first_allele_specific_line = True
first_nonspecific_line = True

for gene_transcript_pos, hits in hit_dict.items():
    found_nas = False
    sums = {'ORE': 0, 'SAM': 0, 'NA': 0}
    for cur_hit in hits:
        sums[cur_hit.background] += fc_count_dict[cur_hit.tag]
        if cur_hit.background == 'NA':
            found_nas = True
    if found_nas == False:
        #Allele specific!
        gene_name = hits[0].gene
        gene_transcript = hits[0].get_id()
        ore_count = str(sums['ORE'])
        sam_count = str(sums['SAM'])
        output_fields = [gene_name, gene_transcript, ore_count, sam_count]
        if first_allele_specific_line:
            first_allele_specific_line = False
        else:
            allelespecific_handle.write("\n")
        allelespecific_handle.write(",".join(output_fields))
    gene_name = hits[0].gene
    gene_transcript = hits[0].get_id()
    count = str(sum(sums.values()))
    output_fields = [gene_name, gene_transcript, count]
    if first_nonspecific_line:
        first_nonspecific_line = False
    else:
        nonspecific_handle.write("\n")
    nonspecific_handle.write(",".join(output_fields))

nonspecific_handle.close()
allelespecific_handle.close()

#The non-allele specific dataset has:
#  gene_name, gene_transcript_start, count
  
#The allele-specific dataset has:
#  gene_name, gene_transcript_start, ORE_count, SAM_count
