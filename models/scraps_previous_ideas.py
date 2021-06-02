from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC  # allows for 'N' in the sequence

from Bio.SeqIO.QualityIO import FastqGeneralIterator


def wt(f, fo):
    n = 40000
    with open(f) as myfile:
        head = [next(myfile) for x in range(n)]

    with open(fo, "w") as fileout:
        fileout.writelines(head)

    seqs = SeqIO.parse(fo, "fastq", IUPAC.ambiguous_dna)
    # print(seqs)

    tmp = []
    for s in seqs:
        tmp.append(Seq(s.seq._data, IUPAC.ambiguous_dna))

    motif = motifs.create(seqs, IUPAC.ambiguous_dna)
    return motif.consensus._data


# partial input file for processing
# compare with WT sequence at this point? Could index
# how fast can I read in a fastq file?

# What do I want to do:
# 1. read all barcodes
# 2. read all variant/wt strings
# 3. associate each barcode with variant/wt string
# 4. count the number of repeated barcode-variant/wt associations
#

# What problems do I have with the current system:
# 1. when I have a new barcode-variant/wt sequencing run, I have to start over and rerun all the scripts
# in a 1-off kind of way. It would be better if I could just add the new sequencing results to the pipeline and
# "turn the crank".
# Should I try to store all reads into a single file?
# Do I need to combine barcode-variant/wt sequencing reads together?
# Maybe I can use each sequencing result separately and just vote for the winner?
# Maybe I can collect all the barcodes into a single file, remove those with less than 10 observations,
# associate the remaining barcodes with variants/wt?



total_len = 0

wt = wt(f_in, f_out)

print(wt)
