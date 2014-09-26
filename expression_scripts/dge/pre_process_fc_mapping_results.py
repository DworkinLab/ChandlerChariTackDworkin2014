'''
Reads a SAM file containing the mappings generated by bwa
Outputs a CSV file that has:
  sequence_tag, gene_name, gene_transcript, start_pos, background(SAM/ORE/NA)
Also outputs another CSV file that has the same information for all the sequence
tags that were problematic during mapping

Usage:
%% python pre_process_fc_mapping_results.py <fc_map.sam> <output.csv> <problems.csv>

'''

import sys

input_handle = open(sys.argv[1], 'r')
output_handle = open(sys.argv[2], 'w')
problem_handle = open(sys.argv[3], 'w')

###############################################
#Class to store information about a mapping hit
###############################################
class mapping_hit:
    def __init__(self, hit_id="", pos=-1, edit_dist=-1):
        self.hit_id = hit_id
        self.pos = pos
        self.edit_dist = edit_dist

#################################################
#Useful functions
#################################################
#Takes a list of strings
#Returns the one that we want, denoted by desired_field
def get_field(fields, desired_field):
    result = None
    for cur_field in fields:
        if cur_field.find(desired_field) == 0:
            result = cur_field.rpartition(":")[2]
    return result

#Takes an XA field
#Returns a list of tuples of the form, each containing (hit_id, position, edit_dist)
def get_hits_from_xa(xa_string):
    xa_string = xa_string.strip() #Remove any trailing whitespace
    xa_string = xa_string.strip(";") #Remove any trailing semicolons
    xa_strings = xa_string.split(";")
    hits = []
    for cur_xa_string in xa_strings:
        cur_hit_fields = cur_xa_string.split(",")
        cur_hit_id = cur_hit_fields[0]
        cur_hit_pos = abs(int(cur_hit_fields[1]))
        cur_hit_edit_dist = int(cur_hit_fields[3])
        cur_hit = mapping_hit(cur_hit_id, cur_hit_pos, cur_hit_edit_dist)
        hits.append(cur_hit)
    return hits

#Takes the best hit reported by bwa as well as all alternative hits
#Returns a list of all hits that are as good as the best one
def get_all_best_hits(best_hit, other_hits):
    best_hits = [best_hit]
    for cur_hit in other_hits:
        #In practice, I've found that BWA may not put the hit with the lowest
        # edit distance as the "best" hit if there are indels
        if cur_hit.edit_dist <= best_hit.edit_dist:
            best_hits.append(cur_hit)
    return best_hits

#A function that allows us to always pick the same transcript/mapping position to report when there are
# reads that map to multiple transcripts of the same gene
def get_lowest_transcript_pos(hits):
    lowest_hit = hits[0]
    for cur_hit in hits:
        if get_transcript(cur_hit.hit_id) < get_transcript(lowest_hit.hit_id):
            lowest_hit = cur_hit
        elif get_transcript(cur_hit.hit_id) == get_transcript(lowest_hit.hit_id):
            if cur_hit.pos < lowest_hit.pos:
                lowest_hit = cur_hit
    return lowest_hit

#Output a hit to a file
def output_hit(seq_tag, hit, handle, background=None):
    output_fields = [seq_tag] #Sequence tag
    output_fields.append(get_gene(hit.hit_id)) #Gene
    output_fields.append(get_transcript(hit.hit_id)) #Transcript ID
    output_fields.append(str(hit.pos)) #Position
    if background == None: #Background -- use the hit's background info if not specified
        output_fields.append(get_background(hit.hit_id))
    else:
        output_fields.append(background)
    output_line = ",".join(output_fields)
    handle.write(output_line + "\n")

#Records an allele-specific hit to the output file
def record_allele_specific_hit(seq_tag, hits):
    #print "Allele specific hit*********************************************"
    #We know all the hits here have the same background
    #So pick the lowest transcript id and go with that
    #Systematically find the lowest transcript/hit
    lowest_hit = get_lowest_transcript_pos(hits)
    output_hit(seq_tag, lowest_hit, output_handle)
    return False

#Records a non-allele-specific hit to the output file
def record_nonspecific_hit(seq_tag, hits):
    #print "Nonspecific hit*************************************************"
    #We know there's just one gene here, and all of the different transcripts (if there are more than one)
    # pair up well
    lowest_hit = get_lowest_transcript_pos(hits)
    output_hit(seq_tag, lowest_hit, output_handle, "NA")
    return False

#Records a hit with problems
def record_problematic_hit(seq_tag, hits):
    #print "Problematic hit*************************************************"
    #Something funny about this sequence tag...
    #It has hits in both SAM and ORE for a single gene, but the locations within that
    # gene that it maps to differs between the two backgrounds
    for cur_hit in hits:
        output_hit(seq_tag, cur_hit, problem_handle)
    return False

#Extracts the gene portion from a hit id
def get_gene(hit_id):
    result = hit_id.partition("-")[2]
    result = result.rpartition("-")[0]
    return result

#Extracts the transcript portion from a hit id
def get_transcript(hit_id):
    result = hit_id.rpartition("-")[2]
    return result

#Extracts the background portion from a hit id
def get_background(hit_id):
    result = hit_id.partition("-")[0]
    return result

#Takes a list of mapping hits
#Goes through them and returns a tuple that contains two lists
#The first list has all the hits that represent complementary SAM-ORE pairs
#The second list has all the hits for which a match could not be found for the other background
def get_matched_and_unmatched_hits(all_hits):
    done = False
    i = 0
    unmatched_hits = all_hits
    matched_hits = []
    while not done:
        cur_hit = unmatched_hits[i]
        #Look for a match to the current hit in the remaining hits
        match_found = False
        for j in range(i+1,len(unmatched_hits)):
            test_hit = unmatched_hits[j]
            matches = (get_gene(cur_hit.hit_id) == get_gene(test_hit.hit_id))
            matches = matches & (get_transcript(cur_hit.hit_id) == get_transcript(test_hit.hit_id))
            matches = matches & (cur_hit.pos == test_hit.pos)
            matches = matches & ( ((get_background(cur_hit.hit_id) == "SAM") & (get_background(test_hit.hit_id) == "ORE")) | ((get_background(cur_hit.hit_id) == "ORE") & (get_background(test_hit.hit_id) == "SAM")) )
            #matches = matches & (((get_background(cur_hit.hit_id) == "SAM") & (get_background(test_hit.hit_id)) == "ORE") | ((get_background(cur_hit.hit_id) == "ORE") & (get_background(test_hit.hit_id) == "SAM")))
            if matches:
                #We have a SAM and ORE match!
                #Pop the matching hits from the list
                #Pop the second one first because popping the first one will screw up the index
                # of the second one
                matched_hits.append(unmatched_hits.pop(j))
                matched_hits.append(unmatched_hits.pop(i))
                match_found = True
                break
        if not match_found:
            i += 1
        if i >= len(unmatched_hits) - 1:
            done = True
    return (matched_hits, unmatched_hits)

##############################################
#The script
##############################################
for cur_line in input_handle:
    #Make sure it's not a header/comment line
    if cur_line[0] != "@":
        fields = cur_line.split("\t")
        cur_tag = fields[0]
        #print cur_tag
        #Make sure there was at least one hit for this read
        if fields[2] != "*":
            #Get the info for the best hit
            best_hit_id = fields[2]
            best_pos = abs(int(fields[3]))
            best_edit_dist = int(get_field(fields, "NM"))
            best_hit = mapping_hit(best_hit_id, best_pos, best_edit_dist)
            #Are there multiple hits?
            xa_field = get_field(fields, "XA")
            if xa_field != None:
                other_hits = get_hits_from_xa(xa_field)
            else:
                other_hits = []
            best_hits = get_all_best_hits(best_hit, other_hits)
            if len(best_hits) == 1:
                #Only one hit, so report it as allele-specific
                record_allele_specific_hit(cur_tag, best_hits)
            else:
                genes = [ get_gene(x.hit_id) for x in best_hits ]
                genes = list(set(genes))
                #Make sure there is only one gene... if there's more than one, we can't do anything with it
                if len(genes) == 1:
                    (matched_hits, unmatched_hits) = get_matched_and_unmatched_hits(best_hits)
                    if len(unmatched_hits) == 0:
                        #There are no unmatched hits -- it must be a non-allele-specific hit
                        record_nonspecific_hit(cur_tag, matched_hits)
                    else:
                        #There ARE unmatched hits
                        #If there are also matched hits, then we can't really have any 
                        #  confidence that this hit is allele-specific, so record it as
                        #  problematic for now
                        #If there are no matched hits, we want to make sure that the hits
                        #  we DO have are all from the same background -- if not, that's a
                        #  problem because it means this hit maps to one location in one
                        #  background and another location in the other
                        if len(matched_hits) != 0:
                            #There are matched hits! This is a problem
                            record_problematic_hit(cur_tag, unmatched_hits)
                            record_problematic_hit(cur_tag, matched_hits)
                        else:
                            #There are unmatched hits but no matched hits -- we want to now
                            #check if the unmatched hits are all from the same background
                            backgrounds = [ get_background(cur_hit.hit_id) for cur_hit in unmatched_hits ]
                            if ("SAM" in backgrounds) & ("ORE" in backgrounds):
                                record_problematic_hit(cur_tag, unmatched_hits)
                            else:
                                record_allele_specific_hit(cur_tag, unmatched_hits)
                        
                        
                        
                    
                        
            

input_handle.close()
output_handle.close()
problem_handle.close()
