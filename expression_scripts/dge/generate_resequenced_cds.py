'''
Uses a flybase CDS sequence (fasta) and re-sequenced fly genomes (fasta) to
make a new CDS sequence incorporating the sequence changes included in 
the re-sequenced genome

Usage:
%% python generate_resequenced_cds.py <flybase_CD.fasta> <resequenced_genome.fa> <output.fasta>

'''

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

#Ignore any genes that aren't in these chromosomes:
valid_chroms = ["X", "2L", "2R", "3L", "3R"]

#Read the original CDS file
cds = []
cds_handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(cds_handle, "fasta"):
    cds.append(record)
cds_handle.close()

#Read the re-sequenced genome
#Comment out during debugging because it's slow
genome_handle = open(sys.argv[2], "rU")
genome = SeqIO.to_dict(SeqIO.parse(genome_handle, "fasta"))
for record in SeqIO.parse(genome_handle, "fasta"):
    genome.append(record)
genome_handle.close()

#Extracts the location field from a sequence description
def extract_loc(description):
    fields = description.split("; ")
    for i in fields:
        if (i[0:4]=="loc="):
            result = i[4:len(i)]
    return(result)

#Returns the modified sequence given the location string
def parse_loc(loc):
    chrom = loc.partition(":")[0]
    if chrom in valid_chroms:
        chunks = loc.partition(":")[2]
        result = parse_chunks(chrom, chunks)
        return result
    else:
        return None
    
#Does all the hard work of parsing the location string
#Note -- need to do it recursively because there may be complement()'s within join()'s!
def parse_chunks(chrom, chunks):
    result = ''
    if (chunks.find("join(") == 0):
        #Currently on a join(...) section so parse the stuff in it
        sub_chunks = chunks.partition("(")[2]
        sub_chunks = sub_chunks.rpartition(")")[0]
        result += parse_chunks(chrom, sub_chunks) 
    elif (chunks.find("complement(") == 0):
        #Currently on a complement(...) section so parse the stuff in it
        sub_chunks = chunks.partition("(")[2]
        sub_chunks = sub_chunks.rpartition(")")[0]
        result += parse_chunks(chrom, sub_chunks).reverse_complement()
    elif ("," in chunks):
        #Currently have a consecutive list of chunks, so parse them
        sub_chunks = chunks.split(",")
        for cur_chunk in sub_chunks:
            result += parse_chunks(chrom, cur_chunk)
    else: #Currently have an individual chunk; parse it
        cur_chunk_fields = chunks.split("..")
        cur_chrom_seq = genome[chrom]
        start = int(cur_chunk_fields[0]) - 1
        stop = int(cur_chunk_fields[1])
        result += cur_chrom_seq[start:stop].seq
    return result

#Gets the percent divergence between two sequences
def get_perc_diff(seq1, seq2):
    if len(seq1) != len(seq2):
        return -1.0
    diff = 0.0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            diff += 1.0
    return diff / len(seq1)

#Gets the number of N's in a sequence (assuming all uppercase)
def get_num_ns(myseq):
    ns = 0
    ns = myseq.count("N")
    return ns

#Replaces all the N's in a sequence with the appropriate base from the template
#Assumes it's already all uppercase
def replace_ns(new_seq, template):
    if len(new_seq) != len(template):
        raise Exception("replace_ns: original and template are not the same length")
    result = new_seq #str(new_seq)
    if "N" in result:
        for i in range(len(result)):
            if result[i] == "N":
                result = result[0:i] + template[i:i+1] + result[i+1:len(result)]
                #result[i] = template[i]
    return result
            

#cds = cds[0:1000]
output_sequences = []
matches = 0
non_matches = 0
count = 0
print "Starting..."
for c in cds:
    if count % 100 == 0:
        print "At number %d" % count
    location = extract_loc(c.description)
    temp_seq = parse_loc(location)
    if temp_seq != None:
        temp_seq = temp_seq.upper() #Make the new sequence all uppercase
        temp_seq = replace_ns(temp_seq, c.seq) #Replace N's in it
        old_string = str(c.seq) #A string version of the original sequence
        new_string = str(temp_seq) #A string version of the new sequence
        new_record = c #A SeqRecord version of the new sequence
        new_record.seq = temp_seq #Replace the sequence in the SeqRecord version of the new sequence
        output_sequences.append(new_record)
        #Check if they are similar enough...
        perc_diff = get_perc_diff(old_string, new_string)
        num_ns = get_num_ns(new_string)
        if (perc_diff > 0.01) | (perc_diff == -1.0):
            non_matches += 1
            print "%d %s %s %s %d %d %f %d" % (count, c.id, old_string[0:10], new_string[0:10], len(old_string), len(new_string), perc_diff, num_ns)
        else:
            matches += 1
    count += 1

output_handle = open(sys.argv[3], 'w')
SeqIO.write(output_sequences, output_handle, "fasta")
output_handle.close()

print "Done! %d sequences matched; %d sequences did not" % (matches, non_matches)