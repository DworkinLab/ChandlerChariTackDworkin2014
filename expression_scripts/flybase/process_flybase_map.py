'''
Takes the flybase map table downloaded directly from flybase and outputs a more user-friendly
version for us.

Output format:
fbgn,gene_symbol,chrom,start,stop

Usage:
%% python process_flybase_map.py <input.tsv> <output.csv>

'''

import sys

input_handle = open(sys.argv[1], 'r')
output_handle = open(sys.argv[2], 'w')

first_line = True
for cur_line in input_handle:
    #Make sure it's not a comment
    cur_line = cur_line.strip()
    if len(cur_line) <= 1:
        continue
    if cur_line[0] == "#":
        continue
    #Make sure it has the right number of fields
    fields = cur_line.split("\t")
    if len(fields) == 5:
        #There really are only four fields, but for some reason one of the tabs
        # is atually a double
        print fields
        gene_symbol = fields[0]
        #Make sure it is actually a D. melanogaster gene
        #If not, its symbol will be in the form "Gspe\gene" where Gspe = Genus species, e.g., "Dsim\Or46a"
        if "\\" in gene_symbol:
            continue
        fbgn = fields[1]
        #Now we need to process the fourth field, which contains all the location information, e.g., 3L:15200944..15201692(1)
        location_string = fields[4]
        if not (":" in location_string):
            continue
        if not (".." in location_string):
            continue
        if not (("(" in location_string) & (")" in location_string)):
            continue
        chrom = location_string.partition(":")[0]
        coords = location_string.partition(":")[2].partition("(")[0]
        start = coords.partition("..")[0]
        stop = coords.partition("..")[2]
        output_fields = [fbgn, gene_symbol, chrom, start, stop]
        if first_line:
            first_line = False
        else:
            output_handle.write("\n")
        output_handle.write(",".join(output_fields))

input_handle.close()
output_handle.close()
