import glob
import os
import models.consensus
import models.labelvariant
import models.c_bcvarfiles
import Bio
import re


# Define the folder where the fastq files are and generate a list of all fastq files
folder = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\"  # input file

file_stem = []
tmp = []
for f in glob.glob(folder + "*fastq"):
    tmp.append(re.split(r"_R",f))
    if re.split(r"_R",f)[0] not in file_stem:
        file_stem.append(re.split(r"_R",f)[0])

suffix = tmp[1][1][1:]

if not os.path.isdir(folder + "data\\"):
    os.mkdir(folder + "data\\")

#  Process raw reads, find 'wt' (or consensus sequence), write out barcode + approximately 10 nucleotides,
#  TODO: find consensus sequence for both reads (R1 and R2), then look for 18 nucleotide stretches where the first
#   10k(?) sequences do not match the consensus.
#  compare each read with 'wt', remove reads with 'N's, write aggregated file to disk
i = 1
for f in file_stem:
    bc_init_R1, r_seq_R1, lobc_R1, robc_R1, bcstart_R1, bcend_R1 = models.consensus.findconsensus(f+"_R1"+suffix)
    bc_init_R2, r_seq_R2, lobc_R2, robc_R2, bcstart_R2, bcend_R2 = models.consensus.findconsensus(f+"_R2"+suffix)
    if bc_init_R1:
        wt = r_seq_R2
        file_barcode = f+"_R1"+suffix
        file_var = f + "_R2" + suffix
        lobc = lobc_R1
        robc = robc_R1
        bcstart = bcstart_R1
        bcend = bcend_R1
    elif bc_init_R2:
        wt = r_seq_R1
        file_barcode = f+"_R2"+suffix
        file_var = f+"_R1"+suffix
        lobc = lobc_R2
        robc = robc_R2
        bcstart = bcstart_R2
        bcend = bcend_R2
    else:
        print("didn't find barcode")
        break
    print(wt+" "+file_barcode+" "+file_var+" "+str(bcend))
    models.consensus.write_bcvar(file_var, file_barcode, folder + "data/" + str(i) + "th_fastq.tmp.wt",
                             folder + "data/" + str(i) + "var", wt, bcstart, bcend, lobc, robc)
    with open(folder + "data/" + str(i) + "th_fastq.wt", 'w') as out_handle:
        out_handle.write(wt)
    i += 1
    break


#  read "wt" and compare with plasmid sequence:
#  check to see if barcode is included in the read? Should this already be removed in the first step?
# I could check for the sequence flanking either side of the barcode and split based on that match?
# I think the step which determines the "wt" sequence should also establish the frame of the sequenced segment.

i = 1
file = folder + "data\\" + str(i) + "var"
result = models.c_bcvarfiles.find_and_reduce_bcvars(file)
result_reads = result.reset_index()

# I want to match the full sequence of KCNH2 with the read and pull out all the variants and calculate HGSVg/p

test = Bio.SeqIO.parse(folder+"data\\plasmid.fasta", 'fasta')
plasmid = next(test)
plasmid = plasmid.upper()
tile_seq = next(test)
tile_seq = tile_seq.upper()
wt_gene = next(test)
wt_gene = wt_gene.upper()
tmp_wt = wt_gene.seq._data
tmp_wt = tmp_wt.upper()
wt_prot = wt_gene.translate()

read_wt = open(folder + "data\\0th_fastq.wt", 'r').readline()
read_wt_seq_for = Bio.Seq.Seq(read_wt)
read_wt_seq_rev = read_wt_seq.reverse_complement()
if tile_seq.seq.find(read_wt_seq_for[0:18]) >= 0 or tile_seq.seq.find(read_wt_seq_for[-18:]) >= 0:
    read_wt = read_wt_seq_for
    reverse = False
elif tile_seq.seq.find(read_wt_seq_rev[0:18]) >= 0 or tile_seq.seq.find(read_wt_seq_rev[-18:]) >= 0:
    read_wt = read_wt_seq_rev
    reverse = True
else:
    print("could not match consensus read to plasmid sequence")



(adj_start, from_beginning, offset) = models.labelvariant.check_outside_cDNA(wt_gene, plasmid, read_wt_seq_rev)
(prot_start, frame_offset) = models.labelvariant.find_frame(read_wt_seq_rev, wt_prot)

result_reads = result_reads.rename(columns={0:'read_count'})

#  TODO: figure out if the sequence needs to be reversed or not. I'll have the 'wt' and
#   the plosmid sequence, so this should be easy...

# TODO: find the differences between single variant and ´wt´, report the variant ID and the location, add the offset
#  and append that character set (´wt´, resnum, and mutAA) to the ´results´ file

if not adj_start:
    for index, row in result_reads.iterrows():
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



