from Bio.Seq import Seq

# List of DNA sequences
sequences = [
    "ATGCGTAC",
    "GCGTTAGC",
    "TTAGGCTA",
    "ATGCTAGG",
    "CGTATGCA"
]

# Function to calculate GC content
def gc_content(seq):
    return round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)

# Function to check sequence validity
def is_valid(seq):
    return all(base in "ATCG" for base in seq)

# Loop through each sequence
for s in sequences:
    print(f"Sequence: {s}")
    print(f"Length: {len(s)}")
    print(f"GC Content: {gc_content(s)}%")
    print(f"Valid: {is_valid(s)}")
    seq_obj = Seq(s)
    print(f"Reverse Complement: {seq_obj.reverse_complement()}")
    print("-" * 30)
