import random
from  collections import Counter
from Bio_structure import Nucleotide_Base, Rna_codons, Dna_codons
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqUtils import GC

class b_seq2:

    def __init__(self, seq="ATCG", seq_type= "DNA", Label= "No Label"):
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.label = Label
        self.isValid = self.__validate()
        assert self.isValid, f"The entered sequence is not a valid {self.seq_type} sequence."

    def __validate(self):
        return set(Nucleotide_Base[self.seq_type]).issuperset(self.seq)

    def gen_rnd_seq(self, length="20", seq_type="DNA"):
        seq = "".join([random.choice(Nucleotide_Base[seq_type]) for i in range(length)])
        self.__init__(seq, seq_type, f"Randomly generated {seq_type} sequence.")
        return self.__init__(seq, seq_type, f"Rnadomly generated {self.seq_type} sequence.")

    def get_info(self):
        return f"ID:\n {self.label}\
            \nSequence Type:\n {self.seq_type}\
            \nLength:\n {len(self.seq)}\
            \nGC-Content:\n {dna.gc()}\
            \nSequence:\n {self.seq}"

    def k_nuc(self, k, f):
        templist =[]
        for i in range(0, len(self.seq)-(k-1), 1):
            templist.append(self.seq[i:i+k])
        frequency = dict(Counter(templist))
        return f" {templist}\n Most Frequent:\n  {Counter.most_common(frequency, f)}"

    def compliment(self):
        seq = Seq(self.seq)
        return seq.complement()

    def reverse_complement(self):
        seq = Seq(self.seq)
        return str(seq.reverse_complement())

    def transcription(self):
        seq = Seq(self.seq)
        return seq.reverse_complement().transcribe()

    def translation(self, init_pos = 0):
        seq = Seq(self.seq)
        return str(seq.translate())

    def gc(self):
        seq = Seq(self.seq)
        return f"{round(GC(seq), 3)}%"

    def nuc_frequency(self):
        return dict(Counter(self.seq))

    def reading_frames(self):
        frames = []
        frames.append(str(self.translation(0)))
        frames.append(str(self.translation(1)))
        frames.append(str(self.translation(2)))
        tempseq = b_seq2(self.reverse_complement(), self.seq_type)
        frames.append(str(tempseq.translation(0)))
        frames.append(str(tempseq.translation(1)))
        frames.append(str(tempseq.translation(2)))
        del(tempseq)
        return frames

    def oriC_hunt(self, k , f):
        templist = []
        for i in range(0, len(self.seq)-(k-1), k):
            templist.append(self.seq[i:i+k])
        frequency = dict(Counter(templist))
        most_freq = Counter.most_common(frequency, f)
        return most_freq

    def scan_rfs(self, aa_seq):
        currentprot = []
        prots = []
        for aa in aa_seq:
            if aa == "_":
                if currentprot:
                    for p in currentprot:
                        prots.append(p)
                    currentprot = []
            else:
                if aa == "M":
                    currentprot.append("")
                for i in range(len(currentprot)):
                    currentprot[i] += aa
        return prots

    def pro_from_rfs(self, startpos = 0, endpos = 0, ordered = False):
        if endpos > startpos:
            tempseq = b_seq2(self.seq[startpos:endpos], self.seq_type)
            rfs = tempseq.reading_frames()
            del tempseq
        else:
            rfs = self.reading_frames()
        res = []
        for rf in rfs:
            prot = self.scan_rfs(str(rf))
            for p in prot:
                res.append(p)
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

for i in SeqIO.parse("sequence.fasta", "fasta"):
    seq = str(i.seq)
    name = i.name

dna = b_seq2(seq, "DNA", name)
#dna.gen_rnd_seq(60, "DNA")


if __name__== "__main__":
    print(dna.get_info())
    print("\nComposition:")
    for i in dna.nuc_frequency():
        print(f" {i} = {dna.nuc_frequency()[i]}")
    print(f"\nDinucleotides:\n {dna.k_nuc(2, 5)}")
    print(f"\nTrinucleotides:\n {dna.k_nuc(3, 5)}")
    print(f"\nHexaNucleotides:\n {dna.k_nuc(6, 5)}")
    print(f"\nNonaNucleotides:\n {dna.k_nuc(9, 5)}")
    print(f"\nCompliment:\n {dna.compliment()}")
    print(f"\nReverse Complement:\n {dna.reverse_complement()}")
    print(f"\nTranscription:\n {dna.transcription()}")
    print(f"\nTranslation:\n {dna.translation()}")
    print("\nOriC hunt:")
    print(f" 2-mers: {dna.oriC_hunt(2,4)}\
        \n 3-mers: {dna.oriC_hunt(3,4)}\
        \n 4-mers: {dna.oriC_hunt(4,4)}\
        \n 6-mers: {dna.oriC_hunt(6,4)}\
        \n 9-mers: {dna.oriC_hunt(9,4)}\n")
    for i in dna.reading_frames():
        print(f"ORF:\n {i}")
    for i in dna.pro_from_rfs():
        print(i)
    for i in dna.reading_frames():
        print(dna.scan_rfs(i))
