'''
Reads the allele-specific combined CSV file
Goes through all the mapping hits
Removes the ones that are probably due to a mutation in the adapter restriction site sequence
Detailed outline:
	-Goes through each "hit" in the original output file (~3800)
	-If it is not SAM- or ORE-specific: keep it (~1200)
	-If only the SAM or ORE allele is expressed, but not both:
		-If no adapter sequence is found in either: discard it (~ 380)
		-If adapter is found in the sequence of the "missing" allele: discard it (~850)
		-If adapter sequence is found in both SAM and ORE sequences: keep it (~1400)

Input/Output format:
gene, fbgn, gene_transcript_start, fc3_ore_count, fc3_sam_count, fc4_ore_count, fc4_sam_count, fc5_ore_count, fc5_sam_count, fc6_ore_count, fc6_sam_count

Usage:
%% python remove_false_allele_specifics.py <allele_specifics.csv> <ore_seq.fasta> <sam_seq.fasta> > output.csv

'''

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

###############################################
#Class to store all our data
###############################################
#Next we need to do the combining...
class gene_transcript_start_record:
    def __init__(self, fbgn="", gene="", gene_transcript_start="", fc3_ore_count=0, fc3_sam_count=0, fc4_ore_count=0, fc4_sam_count=0, fc5_ore_count=0, fc5_sam_count=0, fc6_ore_count=0, fc6_sam_count=0):
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
        
    def get_gene_transcript(self):
        return self.gene_transcript_start.rpartition('-')[0]
    
    def get_start(self):
        return int(self.gene_transcript_start.rpartition('-')[2])

###############################################
# Script body
###############################################

#Reads in the CDS files for ORE and SAM so we can check to see whether the sequences have CATG's
ore_handle = open(sys.argv[2], "rU")
ore_seqs = SeqIO.to_dict(SeqIO.parse(ore_handle, "fasta"))
ore_handle.close()

sam_handle = open(sys.argv[3], "rU")
sam_seqs = SeqIO.to_dict(SeqIO.parse(sam_handle, "fasta"))
sam_handle.close()

num_no_adapter = 0
num_no_adapter_anywhere = 0
tot_seqs = 0
potential_false_hits = 0
adapter_in_both = 0
false_hits = 0

hits_handle = open(sys.argv[1], 'r')
for cur_hit in hits_handle:
    tot_seqs += 1
    cur_fields = cur_hit.strip().split(',')
    #Fields:
    #  gene, fbgn, gene_transcript_start, fc3_ore_count, fc3_sam_count, fc4_ore_count, fc4_sam_count, fc5_ore_count, fc5_sam_count, fc6_ore_count, fc6_sam_count
    #Create a gene_transcript_start_record object for the current line 
    cur_hit_record = gene_transcript_start_record(cur_fields[1], cur_fields[0], cur_fields[2], int(cur_fields[3]), int(cur_fields[4]), int(cur_fields[5]), int(cur_fields[6]), int(cur_fields[7]), int(cur_fields[8]), int(cur_fields[9]), int(cur_fields[10]))
    #Set up some record-keeping...
    real_hit = True
    pres_seq = None
    pres_seq_record = None
    miss_seq_record = None
    #Find out if all of the ORE counts are 0
    if ((cur_hit_record.fc3_ore_count == 0) & (cur_hit_record.fc4_ore_count == 0) & (cur_hit_record.fc5_ore_count == 0) & (cur_hit_record.fc6_ore_count == 0)):
        #This is a SAM-only transcript
        pres_seq_record = sam_seqs[cur_hit_record.get_gene_transcript()]
        miss_seq_record = ore_seqs[cur_hit_record.get_gene_transcript()]
    #Find out if all of the SAM counts are 0
    if ((cur_hit_record.fc3_sam_count == 0) & (cur_hit_record.fc4_sam_count == 0) & (cur_hit_record.fc5_sam_count == 0) & (cur_hit_record.fc6_sam_count == 0)):
        #This is an ORE-only transcript
        pres_seq_record = ore_seqs[cur_hit_record.get_gene_transcript()]
        miss_seq_record = sam_seqs[cur_hit_record.get_gene_transcript()]
    if (pres_seq_record == None):
        #Both SAM and ORE alleles are expressed; keep it
        print cur_hit.strip()
    else:
        #Potentially a false hit; now let's check if it really is
        potential_false_hits += 1
        start_pos = cur_hit_record.get_start()
        #First, check for CATG in it
        pres_seq = str(pres_seq_record.seq[start_pos-1:start_pos+21-1]) #The whole sequence of the transcript that we did detect in the DGE data
        miss_seq = str(miss_seq_record.seq[start_pos-1:start_pos+21-1]) #The whole sequence of the transcript allele that seems to be missing/not expressed
        #If the 'CATG' adapter is in it, then make sure it's near the beginning or end of the read, and check if there's a mutation in the other allele
        if ('CATG' in pres_seq):
            first_found = pres_seq.find('CATG')
            last_found = pres_seq.find('CATG')
            if (first_found == 0):
                #Adapter is at the beginning of the sequence tag
                found_pos = first_found
            elif (last_found >= 16):
                #Adapter is at the end of the sequence tag
                found_pos = last_found
            else:
                #Adapter sequence is present, but not at one of the ends of the sequence
                num_no_adapter += 1
                found_pos = -1
            if (found_pos >= 0): #The adapter is near the beginning or end of the sequence tag
                #Now let's check for a mutation in the missing allele's adapter sequence
                if (miss_seq[found_pos:found_pos+4] == 'CATG'):
                    #Adapter is present in both sequences; keep it!
                    print cur_hit.strip()
                    adapter_in_both += 1
                else:
                    #Adapter sequence is only present in the allele that is expressed
                    #Other allele has a mutation in the adapter and that explains why we didn't detect it
                    false_hits += 1
        else:
            #No adapter sequence in the sequence tag at all! We don't want to keep these, so don't
            # bother printing it
            num_no_adapter_anywhere += 1
            num_no_adapter += 1


#print "%d sequences total" % tot_seqs
#print "%d potential false hits" % potential_false_hits
#print "%d have no adapter at ends" % num_no_adapter
#print "%d have no adapter anywhere" % num_no_adapter_anywhere
#print "%d false hits to discard" % false_hits
#print "%d tags to keep because adapter in both" % adapter_in_both


#Adapter sequence: 
#   CATGNNNNNNNNNNNNNNNNN
#   NNNNNNNNNNNNNNNNNCATG
    
hits_handle.close()