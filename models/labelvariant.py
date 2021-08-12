import math
import numpy as np
import models.c_bcvarfiles
import Bio


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


def find_frame(read_seq, wt_prot):
    n = np.arange(3)
    for j in np.nditer(n):
        read_prot = read_seq[j:].translate()
        if wt_prot.seq.find(read_prot) >= 0:
            return wt_prot.seq.find(read_prot), j.max()
        else:
            continue
    print("Frame for translation not found")


def convert_to_mutation(file_list, folder):
    out_list = []
    for var_file in file_list:
        test = Bio.SeqIO.parse(folder + "plasmid.fasta", 'fasta')
        next(test)
        wt_gene = next(test)
        wt_gene = wt_gene.upper()
        wt_prot = wt_gene.translate()

        read_wt = open(var_file + ".read.wt", 'r').readline()
        read_wt = read_wt.upper()
        read_wt = Bio.Seq.Seq(read_wt)

        (prot_start, frame_offset) = models.labelvariant.find_frame(read_wt, wt_prot)
        result = models.c_bcvarfiles.find_and_reduce_bcvars(var_file)
        result_reads = result.reset_index()
        result_reads = result_reads.rename(columns={0: 'count'})

        count_step = 0
        for index, row in result_reads.iterrows():
            var = row['seq'][frame_offset:]
            read_wt_adj = read_wt[frame_offset:]
            fastq_read = Bio.Seq.Seq(var)

            count = sum(1 for a, b in zip(fastq_read.translate(), read_wt_adj.translate()) if a != b)
            mutAA = "NA"
            resnum = "NA"
            wtAA = "wt"
            var_read_adj = 0
            if count > 1:
                if sum(1 for a, b in zip(fastq_read[1:].translate(), read_wt_adj.translate()) if a != b) == 1:
                    var_read_adj = 1
                    wt_read_adj = 0
                else:
                    variant = "multiple"
            fastq_read_prot = fastq_read[var_read_adj:].translate()
            read_wt_prot = read_wt_adj.translate()

            if count == 1:
                i = 0
                for read, wt in zip(fastq_read_prot._data, read_wt_prot._data):
                    i += 1
                    if read != wt:
                        mutAA = read
                        resnum = i + prot_start
                        wtAA = wt
                        continue
                variant = wtAA + str(resnum) + mutAA
            elif count == 0:
                i = 0
                for read, wt in zip(fastq_read._data, read_wt_adj._data):
                    i += 1
                    if read != wt:
                        i = math.floor(i / 3)
                        mutAA = read_wt_adj[i * 3:i * 3 + 3].translate()._data
                        resnum = i + 1 + prot_start
                        wtAA = mutAA
                        continue
                variant = wtAA + str(resnum) + mutAA
            count_step += 1
            result_reads.at[index, 'variant'] = variant
            result_reads.at[index, 'mutAA'] = mutAA
            result_reads.at[index, 'resnum'] = str(resnum)
            result_reads.at[index, 'wtAA'] = wtAA

        result_reads.to_csv(var_file + 'processed.csv')
        out_list.append(var_file + 'processed.csv')

    return out_list


def count_sorted(file_list):
    for sort_file in file_list:
        result = models.c_bcvarfiles.find_and_reduce_bc_sort(sort_file)
        result_reads = result.reset_index()
        result_reads = result_reads.rename(columns={0: 'count'})
        result_reads.to_csv(sort_file + '.processed.csv')
