import glob
import os
import models.consensus
import models.labelvariant
import models.c_bcvarfiles
import Bio


# Define the folder where the fastq files are and generate a list of all fastq files
folder = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\"  # input file

files = glob.glob(folder + "*fastq")
#files_barcode = glob.glob(folder + "*R2*fastq")
if not os.path.isdir(folder + "data\\"):
    os.mkdir(folder + "data\\")

#  Process raw reads, find 'wt' (or consensus sequence), write out barcode + approximately 10 nucleotides,
#  TODO: find consensus sequence for both reads (R1 and R2), then look for 18 nucleotide stretches where the first
#   10k(?) sequences do not match the consensus.
#  compare each read with 'wt', remove reads with 'N's, write aggregated file to disk
i = 0
for files in files_rev:
    wt = models.consensus.findconsensus(files)
    models.consensus.write_bcvar(files, files_barcode[i], folder + "data/" + str(i) + "th_fastq.tmp.wt",
                                 folder + "data/" + str(i) + "var", wt)
    with open(folder + "data/" + str(i) + "th_fastq.wt", 'w') as out_handle:
        out_handle.write(wt)
    i += 1

#  read "wt" and compare with plasmid sequence:
#  check to see if barcode is included in the read? Should this already be removed in the first step?
# I could check for the sequence flanking either side of the barcode and split based on that match?
# I think the step which determines the "wt" sequence should also establish the frame of the sequenced segment.

i = 0
file = folder + "data/" + str(i) + "th_fastq.tmp"
result = models.c_bcvarfiles.find_and_reduce_bcvars(file)
result_reads = result.reset_index()

# I want to match the full sequence of KCNH2 with the read and pull out all the variants and calculate HGSVg/p

test = Bio.SeqIO.parse(folder+"data\\plasmid.fasta", 'fasta')
plasmid = next(test)
plasmid = plasmid.upper()
wt_gene = next(test)
wt_gene = wt_gene.upper()
tmp_wt = wt_gene.seq._data
tmp_wt = tmp_wt.upper()
wt_prot = wt_gene.translate()

read_wt = open(folder + "data\\0th_fastq.wt", 'r').readline()
read_wt_seq = Bio.Seq.Seq(read_wt)
read_wt_seq_rev = read_wt_seq.reverse_complement()

(adj_start, from_beginning, offset) = models.labelvariant.check_outside_cDNA(wt_gene, plasmid, read_wt_seq_rev)
(prot_start, frame_offset) = models.labelvariant.find_frame(read_wt_seq_rev, wt_prot)


result_reads_var_ind = result_reads['wt'] != 'wt'
result_reads_var = result_reads[result_reads_var_ind]
result_reads_var = result_reads_var.rename(columns={0:'read_count'})

# TODO: find the differences between single variant and ´wt´, report the variant ID and the location, add the offset
#  and append that character set (´wt´, resnum, and mutAA) to the ´results´ file

if not adj_start:
    for index, row in result_reads_var.iterrows():
        bc = row['bc']
        read_count = row['read_count']
        var = row['seq']
        fastq_read = Bio.Seq.Seq(var)
        fastq_read_r = fastq_read.reverse_complement()
        fastq_read_r = fastq_read_r[frame_offset:]
        print("bc")
        fastq_read_r_prot = fastq_read_r.translate()

        count = sum(1 for a, b in zip(fastq_read_r_prot, read_wt_seq_rev.translate()) if a != b)
        #break
        mutAA = "NA"
        resnum = "NA"
        wtAA = "wt"
        if count > 1:
            variant = "multiple"
        elif count == 1:
            variant = "Good!"
            i = 0
            for read, wt in zip(fastq_read_r_prot._data,read_wt_seq_rev.translate()._data):
                i += 1
                if read != wt:
                    print(read)
                    mutAA = read
                    resnum = i + prot_start
                    wtAA = wt
                    continue
        else:
            variant = "wt"
        print(fastq_read_r_prot+" "+variant+" "+wtAA+str(resnum)+mutAA)


    if from_beginning:

    elif not from_beginning:
        stuff
    else:
        print("Something informative")

        #  take each read, reverse complement it if needed, then start at nucleotide <offset> and translate the coding
#  sequence, compare to true WT, note the variant and associate the variant protein coding, cDNA variant,
#  and the codon with the barcode



