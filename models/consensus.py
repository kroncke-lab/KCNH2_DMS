from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC  # allows for 'N' in the sequence
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def findconsensus(file_in):
    tmp = []
    with open(file_in) as in_handle:
        count = 0
        for title, seq, qual in FastqGeneralIterator(in_handle):
            tmp.append(Seq(seq, IUPAC.ambiguous_dna))
            count += 1
            if count == 10000:
                wt = motifs.create(tmp, IUPAC.ambiguous_dna).consensus
                break
    return wt._data


def write_bcvar(file_in, file_barcode, file_var_out, file_wt_out, wildtype):
    with open(file_in, 'r') as f, open(file_var_out, 'w') as fout_var, open(file_wt_out, 'w') as fout_wt, open(
            file_barcode, 'r') as barcode_handle:
        fastqline = 0
        for (f_id, seq, q), (r_id, bc_seq, r_q) in zip(FastqGeneralIterator(f), FastqGeneralIterator(barcode_handle)):
            bc = bc_seq[20:54]  # collect a couple nucleotides before and after barcode
            fastqline += 1
            wt_bool = "variant"
            if any(e in seq for e in 'N'):  # remove reads with 'N's
                continue
            if any(e in bc_seq for e in 'N'):  # remove reads with 'N's
                continue
            if seq == wildtype:
                wt_bool = "wt"
                fout_wt.write("{barcode},{seq}\n".format(barcode=bc, seq=seq))
                continue
            fout_var.write("{barcode},{seq}\n".format(barcode=bc, seq=seq))
            if fastqline / 100000.0 == round(fastqline / 100000.0):
                print(fastqline)
            if fastqline == 5000000:
                print("Well, we're close!")
