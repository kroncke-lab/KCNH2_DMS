import os
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC  # allows for 'N' in the sequence
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import glob


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


def write_bcvar(file_in, file_barcode, file_out, wildtype):
    with open(file_in, 'r') as f, open(file_out, 'w') as fout, open(file_barcode, 'r') as barcode_handle:
        fastqline = 0
        for (f_id, seq, q), (r_id, bc_seq, r_q) in zip(FastqGeneralIterator(f), FastqGeneralIterator(barcode_handle)):
            bc = bc_seq[20:54]  # collect a couple nucleotides before and after barcode
            fastqline += 1
            wt_bool = "variant"
            if any(e in seq for e in 'N'):  # remove reads with 'N's
                continue
            if seq == wildtype:
                wt_bool = "wt"
            if fastqline / 100000.0 == round(fastqline / 100000.0):
                print(fastqline)
            fout.write("{barcode},{wt_bool},{seq}\n".format(barcode=bc, wt_bool=wt_bool, seq=seq))
            if fastqline == 5000000:
                print("Well, we're close!")


# Define the folder where the fastq files are and generate a list of all fastq files
folder = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\"  # input file

files_rev = glob.glob(folder + "*R1*fastq")
files_barcode = glob.glob(folder + "*R2*fastq")
if not os.path.isdir(folder+"data\\"):
    os.mkdir(folder+"data\\")

i = 0
for files in files_rev:
    wt = findconsensus(files)
    write_bcvar(files, files_barcode[i], folder+"data/"+str(i)+"th_fastq.tmp", wt)
    with open(folder+"data/"+str(i)+"th_fastq.wt", 'w') as out_handle:
        out_handle.write(wt)
    i += 1


def reduce_fastq(chunk):
    df_min = chunk.groupby(chunk.columns.tolist(), as_index=False).size()
    df_filt = df_min[df_min['size'] > 1]
    df_filt_nonwt = df_filt[df_filt['wt'] != 'wt']
    return df_filt_nonwt

# Count unique barcode-variant pairs, discard n=1 observations (maybe also remove n=2, too?), save file 'to_feather'?
i = 0
file = folder + "data/" + str(i) + "th_fastq.tmp"
reduce_fastq(file)

df_chunks = pd.read_csv(file, names=["bc", "wt", "seq"], chunksize=2000000)

    frames = [reduce_fastq(chunk) for chunk in df_chunks]
pd.df_filt_nonwt.to_feather

