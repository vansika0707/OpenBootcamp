# find_orfs.py
from Bio.Seq import Seq

def get_orfs(sequence):
    sequence = Seq(sequence)
    orfs = set()

    # Forward frames
    for frame in range(3):
        trans = sequence[frame:].translate(to_stop=False)
        for prot in trans.split("*"):
            if prot.startswith("M") and "M" in prot:
                orfs.add(prot)

    # Reverse frames
    rev_seq = sequence.reverse_complement()
    for frame in range(3):
        trans = rev_seq[frame:].translate(to_stop=False)
        for prot in trans.split("*"):
            if prot.startswith("M") and "M" in prot:
                orfs.add(prot)

    return orfs


# Read the DNA sequence, skipping FASTA headers (lines starting with '>')
with open("test_dna_orf.txt") as f:
    lines = f.readlines()
    dna_seq = "".join(line.strip() for line in lines if not line.startswith(">"))

#Get unique ORFs
unique_orfs = get_orfs(dna_seq)

# Print the result
print(f"Found {len(unique_orfs)} unique ORFs:\n")
for i, orf in enumerate(unique_orfs, 1):
    print(f">ORF_{i}\n{orf}")
