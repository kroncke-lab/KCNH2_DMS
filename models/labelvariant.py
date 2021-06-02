import Bio.Seq
from Bio import SeqIO

# I want to match the full sequence of KCNH2 with the read and pull out all the variants and calculate HGSVg/p
folder = 'E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\'

test = SeqIO.parse(folder+"data\\plasmid.fasta", 'fasta')
plasmid = next(test)
wt_gene = next(test)
tmp_wt = wt_gene.seq._data
tmp_wt = tmp_wt.upper()
wt_prot = wt_gene.translate()

read = open(folder + "data\\2th_fastq.wt", 'r').readline()

read_seq = Bio.Seq.Seq(read)
read_seq_rev = read_seq.reverse_complement()

