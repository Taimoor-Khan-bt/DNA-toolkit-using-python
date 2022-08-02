from Bio_sequence import bio_seq
from Bio_structure import Dna_codons, covid_spike
from Bio import SeqIO

#file = SeqIO.parse("sequence.fasta", "fasta")
for i in SeqIO.parse("cov_sequence.fasta", "fasta"):
    file1 = str(i.seq)
    id = i.id

test_dna2 = bio_seq(file1, "DNA", id)

#print(test_dna.gen_random_seq(40, "DNA"))
print(test_dna2.get_seq_info())
print("Nucleotide Frequency:")
for i in test_dna2.nuc_frequency():
    print(f"{i} = {test_dna2.nuc_frequency()[i]}")
print(f"GC-Content: {test_dna2.gc_content()}")
#print(f"Subsequent GC_Content: {test_dna2.subsec_gc_content()}")
print(f"\nDinucleotides: \n{test_dna2.di_nuc()}")
print(f"\nTrinucleotides: \n{test_dna2.tri_nuc()}")
print(f"\nHexanucleotides: \n{test_dna2.hexa_nuc()}")
#print(f"\nTranscription: \n{test_dna2.transcription()}")
print(f"\nReverse Compliment: \n{test_dna2.reverse_compliment()}")
#print(f"\nTranslation: \n{test_dna2.translate()}")
#print(test_dna.codon_usage("L"))
n = 0
for frame in test_dna2.gen_reading_frames():
    n += 1
    print(f"\nORF-{n}: {frame}")
print()
n = 0

#print(f"Scan ORFs: {test_dna2.scan_rf()}")
#print(test_dna2.proteins_from_orfs(ordered=True))
for i in test_dna2.proteins_from_orfs(ordered=True):
    n += 1
    print(f"Possible Protein[{n}][{len(i)}]: {i}")

