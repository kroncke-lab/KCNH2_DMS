import numpy as np


# here I should compare the plasmid seq. with the read 'wt' if DNA in 'wt' doesn't match the gene
def check_outside_cDNA(wt_gene, plasmid, test_seq):
    adj_start = False
    from_beginning = True
    wt_gene = wt_gene.upper()
    plasmid = plasmid.upper()
    test_seq = test_seq.upper()
    if wt_gene.seq.find(test_seq) < 0:
        tmp = plasmid.seq.split(test_seq)
        adj_start = True
        #  the 6,000 below comes from the place of KCNH2 in the full plasmid: there are 3,118 nucleotides before KCNH2
        #  starts, so 6,000 is safe since here I only need to tell if its the sequencing goes beyond the beginning or
        #  follows after the end
        if (len(tmp[0]._data) < 6000):
            from_beginning = True
            seq_use = wt_gene.seq.split(tmp[1]._data[0:200])[0]
            len_use = len(test_seq.split(seq_use)[0]._data)
            return adj_start, from_beginning, len_use
        elif (len(tmp[0]._data) > 6000):
            from_beginning = False
            seq_use = wt_gene.seq.split(tmp[0]._data[-200:])[1]
            len_use = len(test_seq.split(seq_use)[1]._data)
            return adj_start, from_beginning, len_use
        else:
            raise ValueError(
                "Sequence given not found in plasmid provided."
            )
    else:
        return adj_start, from_beginning, 0

def find_frame(read_seq_rev, wt_prot):
    n = np.arange(3)
    for j in np.nditer(n):
        read_prot = read_seq_rev[j:].translate()
        if (wt_prot.seq.find(read_prot) > 0):
            return wt_prot.seq.find(read_prot), j.max()
        else:
            continue
    print("Frame for translation not found")

