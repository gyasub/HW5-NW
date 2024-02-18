# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    # Creating an instance of the NW class
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', float(-10), float(-1))

    # Creating empty list for score ranking

    # Creating empty list for score ranking
    rank = []

    # Aligning each species BRD2 sequence to the human BRD2 sequence
    species_seqs = {
        'Gallus gallus': gg_seq,
        'Mus musculus': mm_seq,
        'Balaeniceps rex': br_seq,
        'Tursiops truncatus': tt_seq
    }

    for species, seq in species_seqs.items():
        output = nw.align(hs_seq, seq)
        score, _, _ = output
        rank.append((score, species))

    # Sorting the species by alignment score in descending order
    rank.sort(reverse=True)

    # Printing species in order of most similar to human BRD2
    print("Species in order of most similar to human BRD2:")
    for score, species in rank:
        print(f"{species}: Alignment score = {score}")
    
    
    # Printing all alignment scores between each species BRD2 and human BRD2
    print("\nAlignment scores between each species BRD2 and human BRD2:")
    for species, seq in species_seqs.items():
        output = nw.align(hs_seq, seq)
        score, _, _ = output
        print(f"{species}: Alignment score = {score}")

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # Printing all alignment scores between each species BRD2 and human BRD2
    
    print("\nAlignment scores between each species BRD2 and human BRD2:")

    for species, seq in species_seqs.items():
        output = nw.align(hs_seq, seq)
        score, _, _ = output
        print(f"{species}: Alignment score = {score}")
    
if __name__ == "__main__":
    main()

