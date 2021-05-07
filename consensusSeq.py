from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC # allows for 'N' in the sequence

f = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\5406-RU-4_S01_L005_R1_001.fastq"  # input file
fo = 'C:\\Users\\KRONCKE\\Box Sync\\Kroncke_lab\\kcnh2\\sequencing\\VANTAGE\\KCNH2_DMS_project\\data\\tmp.fastq'  #
# partial input file for processing

N = 40000
with open(f) as myfile:
    head = [next(myfile) for x in range(N)]


with open(fo, "w") as fileout:
    fileout.writelines(head)

seqs = SeqIO.(fo, "fastq",  IUPAC.ambiguous_dna)
#print(seqs)

tmp = []
for s in seqs:
    tmp.append(Seq(s.seq._data, IUPAC.ambiguous_dna))

motif = motifs.create(seqs, IUPAC.ambiguous_dna)
print(motif.consensus._data)
