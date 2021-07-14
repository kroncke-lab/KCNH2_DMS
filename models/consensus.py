from Bio.Seq import Seq
from Bio import motifs
import re
from Bio.Alphabet import IUPAC  # allows for 'N' in the sequence
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd


def check_align(sequence, data, index):
    unmatched = 0
    for line in data:
        if not re.search(sequence[index:index + 18], line):
            unmatched = unmatched + 1
    return unmatched


def findconsensus(file_in):
    tmp = []
    data = []
    appended_data = []
    prev_match = 0
    diff = 0
    cumulative = 0
    barcode = False
    lobc = None
    robc = None
    bcstart = None
    bcend = None
    with open(file_in) as in_handle:
        count = 0
        for title, seq, qual in FastqGeneralIterator(in_handle):
            tmp.append(Seq(seq, IUPAC.ambiguous_dna))
            data.append(seq)
            count += 1
            if count == 40000:
                wt = motifs.create(tmp, IUPAC.ambiguous_dna).consensus
                break

    pbc = wt._data

    for ind in range(1, 60, 1):
        match = check_align(pbc, data, ind)
        if prev_match > 0:
            diff = match - prev_match
        prev_match = match
        row = [ind, match, diff]
        appended_data.append(row)
        cumulative = cumulative + match
        if cumulative > 720000:
            bc = pbc
            barcode = True
            df = pd.DataFrame(appended_data, columns=["index", "matches", "diff"])
            upper = df[['diff']].idxmin()
            lobc = bc[upper.values[0] - 18 - 5:upper.values[0] - 18 - 1]
            robc = bc[upper.values[0] + 1:upper.values[0] + 5]
            bcstart = upper.values[0]-18-7
            bcend = upper.values[0]+7

    return barcode, wt._data, lobc, robc, bcstart, bcend


def write_bcvar(file_in, file_barcode, file_var_out, file_wt_out, wildtype, bc_start, bc_end):
    with open(file_in, 'r') as f, open(file_var_out, 'w') as fout_var, open(file_wt_out, 'w') as fout_wt, open(
            file_barcode, 'r') as barcode_handle:
        fastqline = 0
        for (f_id, seq, q), (r_id, bc_seq, r_q) in zip(FastqGeneralIterator(f), FastqGeneralIterator(barcode_handle)):
            bc = bc_seq[bc_start:bc_end]  # collect a couple nucleotides before and after barcode
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
