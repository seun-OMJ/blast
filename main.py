from Bio.Align import substitution_matrices
from itertools import product
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import AhoCorasick as AC
blosum62 = substitution_matrices.load("BLOSUM62")

fasta_file = sys.argv[1]
q_seq = sys.argv[2]
k_mer = int(sys.argv[3])
tScore = int(sys.argv[4])
eT = int(sys.argv[5])
max_s = int(sys.argv[6])

kmers = []
neighbourhood = []
kmer_neighbourhood = {}
start_time = time.time()
# 1.3 PRODUCE AND PRINT THE KMERS
if 0 < k_mer < len(q_seq):
    for i in range(len(q_seq) - k_mer + 1):
        subsequence = q_seq[i:i + k_mer]
        kmers.append(subsequence)
print("kmers:", kmers)



def calculate_similarity_score(seq1, seq2):
    blosum_score = 0
    for res1, res2 in zip(seq1, seq2):
        blosum_score += blosum62[res1][res2]
    return blosum_score

sequences = [record for record in SeqIO.parse(fasta_file, "fasta")]




#1.4 CALCULATE AND PRINT THE NEIGHBOURHOOD
alphabet = "ABCDEFGHIKLMNPQRSTVWWXYZ"
three_letter_permutations = product(alphabet, repeat=k_mer)
permutation_strings = [''.join(perm) for perm in three_letter_permutations]

for kmer in kmers:
    for perm in permutation_strings:
        score = calculate_similarity_score(kmer, perm)
        if score >= tScore:
            kmer_neighbourhood[perm] = kmer
            neighbourhood.append(perm)
hit_pos = {}

print("neighbourhood:", neighbourhood)



# 1.5: SEARCHING FOR THE SEEDS
ac_preprocess = AC.preprocess(neighbourhood)
for seq in sequences:
    hit = AC.search(ac_preprocess, str(seq.seq))
    hit_pos[seq.description] = hit





# 1.6: EXTEND THE SEED
def extend_seed(seed,  database_seq, query, eT,start_index):
    key = kmer_neighbourhood.get(seed,"")    
    
    max_score = calculate_similarity_score(key, seed)
    extended_seed = seed

    
    db_indices = -1
    key_index = query.index(key)
    
    # extend to the right
    for i in range(key_index+k_mer , len(query)):
        extended_seed = extended_seed + query[i]
        aligned_seq = database_seq[start_index: start_index +len(extended_seed)]
        score = calculate_similarity_score(extended_seed, aligned_seq)
        
        if score > max_score:
            max_score = score

        if max_score - score > eT:
            extended_seed = extended_seed[:-1]
            break

        

    # extend to the left
    query_indices = -1
    preext = len(extended_seed)
    for i in range(key_index, -1, -1):
        extended_seed = query[i] + extended_seed
        aligned_seq = database_seq[(start_index - key_index + i):start_index+preext]
        
        score = calculate_similarity_score(extended_seed, aligned_seq)
        
        if score > max_score:
            max_score = score
            
        if max_score - score > eT:
            extended_seed = extended_seed[1:-1]
            query_indices = i
            db_indices = start_index - key_index + i
            break

        

    return extended_seed, max_score, db_indices, query_indices






# 1.7: REPORT HSPs WITH A SCORE >=S
hsp_count = 0
next_seq = False
for sequence in sequences:
    next_seq = False
    for seed in hit_pos.get(sequence.description, []):
        if next_seq == True:
            break
        loop_hole_for_seed_dict = hit_pos.get(sequence.description, [])
        for key,value in loop_hole_for_seed_dict.items():
            if next_seq == True:
                break
            for i in range(len(value)):
                result = extend_seed(key, str(sequence.seq), q_seq, eT, value[i])
                
                if result[1] >= max_s:
                    extended_seed, max_score, db_indices, query_indices = result
                    next_seq = True
                    hsp_count += 1
                    print(f"Found one match with score {max_score} in sequence: {sequence.description}")                    
                    print(f"{sequence.seq[db_indices:db_indices + len(extended_seed)]} from the sequence between positions[{db_indices}, {db_indices + len(extended_seed)-1}]")
                    print(f"{extended_seed} from the sequence between positions[{query_indices}, {query_indices + len(extended_seed)-1}]")
                    print()
                    i = len(value) + 1             
end_time = time.time()
execution_time = end_time - start_time
print(f"In total, {hsp_count} hits were found in the database")
print(f"The search was completed in : {execution_time} seconds")

