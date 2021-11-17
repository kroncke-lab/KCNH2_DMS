#  TODO: check that the barcode sequence is reversed when necessary (only necessary for subassembly)



from Bio.Seq import Seq
from Bio import motifs
import Bio
import re
from Bio.Alphabet import IUPAC  # allows for 'N' in the sequence
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd


def check_align(sequence, data, index, step_size):
    unmatched = 0
    for line in data:
        if not re.search(sequence[index:index + step_size], line):
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
            if not re.search('GGGGGGGGGGGGGGGGG', seq):
                data.append(seq)
                count += 1
            if count == 40000:
                wt = motifs.create(tmp, IUPAC.ambiguous_dna).consensus
                break

    pbc = wt._data

    for ind in range(1, 60, 1):
        match = check_align(pbc, data, ind, 18)
        if prev_match > 0:
            diff = match - prev_match
        prev_match = match
        row = [ind, match, diff]
        appended_data.append(row)
        cumulative = cumulative + match
    m = max(appended_data, key=lambda x: x[1])
    m_max = m[1]

    if cumulative > 720000 and m_max > 35000:
        bc = pbc
        barcode = True
        df = pd.DataFrame(appended_data, columns=["index", "matches", "diff"])

        upper = df[['diff']].idxmin()
        bcstart = upper.values[0]-18-9
        bcend = upper.values[0]+10

        appended_data = []
        for ind in range(bcstart, bcend, 1):
            unmatched = 0
            for line in data:
                if not re.search(bc[ind], line[ind]):
                    unmatched = unmatched + 1
            row = [ind, unmatched/len(data)]
            appended_data.append(row)
        for f in appended_data:
            if f[1] > 0.4:
                lobc_end = f[0]
                break
        lobc = bc[lobc_end - 6:lobc_end]
        robc = bc[lobc_end + 18:lobc_end + 18 + 6]

    return barcode, wt._data, lobc, robc, bcstart, bcend


def write_bcvar(file_in, file_barcode, file_var_out, file_wt_out, wildtype, bc_start, bc_end, lobc, robc, var_start,
                var_end, reverse):

    with open(file_in, 'r') as f, open(file_var_out, 'w') as fout_var, open(file_wt_out, 'w') as fout_wt, open(
            file_barcode, 'r') as barcode_handle:

        fastqline = 0

        for (f_id, seq, q), (r_id, bc_seq, r_q) in zip(FastqGeneralIterator(f), FastqGeneralIterator(barcode_handle)):
            fastqline += 1
            bc = bc_seq[bc_start:bc_end]  # collect a couple nucleotides before and after barcode
            m = re.split(lobc, bc[0:11])
            bc = bc[len(m[0]) + len(lobc):len(m[0]) + len(lobc) + 18]

            seq = Bio.Seq.Seq(seq)
            bc = Bio.Seq.Seq(bc)
            if reverse:
                seq = seq.reverse_complement()
            else:
                bc = bc.reverse_complement()
            seq = seq[var_start:var_end]

            if any(e in seq for e in 'N'):  # remove reads with 'N's
                continue

            if seq == wildtype:
                fout_wt.write("{barcode},{seq}\n".format(barcode=bc, seq=seq))
                continue

            fout_var.write("{barcode},{seq}\n".format(barcode=bc, seq=seq))

            if fastqline / 100000.0 == round(fastqline / 100000.0):
                print(fastqline)

            if fastqline == 5000000:
                print("Well, we're close!")


def write_bc_sorted_cells(file_barcode, file_bc_out, bc_start, bc_end, lobc, robc):

    with open(file_bc_out, 'w') as fout_bc, open(file_barcode, 'r') as barcode_handle:

        fastqline = 0

        for (r_id, bc_seq, r_q) in FastqGeneralIterator(barcode_handle):
            fastqline += 1
            bc = bc_seq[bc_start:bc_end]  # collect a couple nucleotides before and after barcode
            m = re.split(lobc, bc[0:11])
            bc = bc[len(m[0]) + len(lobc):len(m[0]) + len(lobc) + 18]

            if any(e in bc for e in 'N'):  # remove reads with 'N's
                continue

            fout_bc.write("{barcode}\n".format(barcode=bc))

            if fastqline / 100000.0 == round(fastqline / 100000.0):
                print(fastqline)

            if fastqline == 5000000:
                print("Well, we're close!")
