import glob
import os
import models.consensus

# Define the folder where the fastq files are and generate a list of all fastq files
folder = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\"  # input file

files_rev = glob.glob(folder + "*R1*fastq")
files_barcode = glob.glob(folder + "*R2*fastq")
if not os.path.isdir(folder + "data\\"):
    os.mkdir(folder + "data\\")

#  Process raw reads, find 'wt' (or consensus sequence), write out barcode + approximately 10 nucleotides,
#  compare each read with 'wt', remove reads with 'N's, write aggregated file to disk
i = 0
for files in files_rev:
    wt = models.consensus.findconsensus(files)
    models.consensus.write_bcvar(files, files_barcode[i], folder + "data/" + str(i) + "th_fastq.tmp.wt",
                                 "data\\var" + str(i), wt)
    with open(folder + "data/" + str(i) + "th_fastq.wt", 'w') as out_handle:
        out_handle.write(wt)
    i += 1
