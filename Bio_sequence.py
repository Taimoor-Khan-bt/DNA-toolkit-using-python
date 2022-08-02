from Bio_structure import Nucleotide_Base, Dna_codons, Rna_codons
import random
from collections import Counter

class bio_seq:

    def __init__(self, seq = "ATCG", seq_type = "DNA", label = "No Label"):
        self.seq = seq
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    def __validate(self):
        """Validates the given DNA sequence"""
        return set(Nucleotide_Base[self.seq_type]).issuperset(self.seq)

    def get_seq_type(self):
        """It gets sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """It gets sequence information"""
        return f"Label : {self.label}\nSequence : {self.seq}\nSequence Type : {self.seq_type}\nLength : {len(self.seq)}"

    def gen_random_seq(self, length= 20, seq_type = "DNA"):
        """Generates random DNA sequence"""
        seq = "".join([random.choice(Nucleotide_Base[seq_type]) for i in range(length)])
        self.__init__(seq, seq_type, f"Randomly generated {seq_type} sequecence")

    def nuc_frequency(self):
        return dict(Counter(self.seq))

    def transcription(self):
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence!"

    def gc_content(self):
        return round((self.seq.count('G') + self.seq.count('C'))/ len(self.seq), 3)

    def reverse_compliment(self):
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def subsec_gc_content(self, k = 20):
        res = []
        for i in range(0, len(self.seq)-k+1, k):
            subsec = self.seq[i:i+k]
            res.append(round((subsec.count('G') +
                             subsec.count('C') )/ len(subsec), 3))
        return res

    def translate(self, init_pos = 0):
        if self.seq_type == "DNA":
            return [Dna_codons[self.seq[i:i + 3]] for i in range(init_pos, len(self.seq)-2, 3)]
        elif self.seq_type == "RNA":
            return [Rna_codons[self.seq[i:i + 3]] for i in range(init_pos, len(self.seq)-2, 3)]

    def codon_usage(self, aminoacid):
        """Provides Frequency of the codon encoding a given aminoacid in a DNA sequence"""
        tempList = []

        if self.seq_type == "DNA":
            for i in range(0, len(self.seq)-2, 3):
                if Dna_codons[self.seq[i:i + 3]] == aminoacid:
                    tempList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq)-2, 3):
                if Rna_codons[self.seq[i:i + 3]] == aminoacid:
                    tempList.append(self.seq[i:i + 3])

        freqDic = dict(Counter(tempList))
        totalWeight = sum(freqDic.values())
        for seq in freqDic:
            freqDic[seq] = round(freqDic[seq] / totalWeight, 2)
        return freqDic

    def gen_reading_frames(self):
        frames = []
        frames.append(self.translate(0))
        frames.append(self.translate(1))
        frames.append(self.translate(2))
        temp_seq = bio_seq(self.reverse_compliment(), self.seq_type)
        frames.append(temp_seq.translate(0))
        frames.append(temp_seq.translate(1))
        frames.append(temp_seq.translate(2))
        del temp_seq
        return frames

    def scan_rf(self, aa_seq):
        """Scan the reading frames for potential proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                #STOPS accumulating aa if stop codon is detected
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    #  START accumulating protein when start codon is detected
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def di_nuc(self):
        tempseq = []
        for i in range(0, len(self.seq)-1, 1):
            tempseq.append(self.seq[i:i+2])
        return f"Length: {len(tempseq)}\n {tempseq}"

    def tri_nuc(self):
        tempseq = []
        for i in range(0, len(self.seq)-2, 1):
            tempseq.append(self.seq[i:i+3])
        return f"Length: {len(tempseq)}\n {tempseq}"

    def hexa_nuc(self):
        tempseq = []
        for i in range(0, len(self.seq)-5, 1):
            tempseq.append(self.seq[i:i+6])
        return f"Length: {len(tempseq)}\n {tempseq}"

    def proteins_from_orfs(self, startpos=0, endpos=0, ordered=False):

        if endpos > startpos:
            temp_seq = bio_seq(self.seq[startpos: endpos], self.seq_type)
            rfs = temp_seq.gen_reading_frames()
            del temp_seq
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.scan_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
                return sorted(res, key=len, reverse=True)
        return res

    
