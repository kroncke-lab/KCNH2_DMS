import glob
import re
import models.consensus
import Bio


def file_grab(folder):
    file_stem = []
    tmp = []
    for f in glob.glob(folder + "*fastq"):
        tmp.append(re.split(r"_R", f))
        if re.split(r"_R", f)[0] not in file_stem:
            file_stem.append(re.split(r"_R", f)[0])

    suffix = tmp[1][1][1:]

    return suffix, file_stem


def file_clean(folder, file_stem, suffix):

    var_file_list = []

    for f in file_stem:
        print(f)
        bc_init_R1, r_seq_R1, lobc_R1, robc_R1, bcstart_R1, bcend_R1 = models.consensus.findconsensus(
            f + "_R1" + suffix)
        bc_init_R2, r_seq_R2, lobc_R2, robc_R2, bcstart_R2, bcend_R2 = models.consensus.findconsensus(
            f + "_R2" + suffix)

        if (bc_init_R1 and not bc_init_R2) or (bc_init_R2 and bc_init_R1 and bcstart_R1 < bcstart_R2):
            wt = r_seq_R2
            file_barcode = f + "_R1" + suffix
            file_var = f + "_R2" + suffix
            lobc = lobc_R1
            robc = robc_R1
            bcstart = bcstart_R1
            bcend = bcend_R1
        elif (bc_init_R2 and not bc_init_R1) or (bc_init_R2 and bc_init_R1 and bcstart_R2 < bcstart_R1):
            wt = r_seq_R1
            file_barcode = f + "_R2" + suffix
            file_var = f + "_R1" + suffix
            lobc = lobc_R2
            robc = robc_R2
            bcstart = bcstart_R2
            bcend = bcend_R2
        else:
            print("Problem finding barcode")
            continue

        var_file = f + "var"
        var_file_list.append(var_file)
        wt_file = f + "th_fastq.tmp.wt"
        wt_seq_file = var_file + ".read.wt"

        wt = Bio.Seq.Seq(wt)
        var_start, var_end, reverse = models.file_structure.clip_read(folder, wt)

        if reverse:
            wt = wt.reverse_complement()
        wt = wt[var_start:var_end]

        with open(wt_seq_file, 'w') as out_handle:
            out_handle.write(str(wt))

        models.consensus.write_bcvar(file_var, file_barcode, var_file, wt_file, wt, bcstart, bcend, lobc, robc,
                                     var_start, var_end, reverse)

    return var_file_list


def sorted_file_clean(file_stem, suffix):

    bc_file_list = []

    for f in file_stem:
        print(f)
        bc_init_R1, r_seq_R1, lobc_R1, robc_R1, bcstart_R1, bcend_R1 = models.consensus.findconsensus(
            f + "_R1" + suffix)
        bc_init_R2, r_seq_R2, lobc_R2, robc_R2, bcstart_R2, bcend_R2 = models.consensus.findconsensus(
            f + "_R2" + suffix)

        if bc_init_R1:
            file_barcode = f + "_R1" + suffix
            lobc = lobc_R1
            robc = robc_R1
            bcstart = bcstart_R1
            bcend = bcend_R1
        elif bc_init_R2:
            file_barcode = f + "_R2" + suffix
            lobc = lobc_R2
            robc = robc_R2
            bcstart = bcstart_R2
            bcend = bcend_R2
        else:
            print("didn't find barcode")
            continue

        file_bc_out = f + ".sorted_bc"
        bc_file_list.append(file_bc_out)

        models.consensus.write_bc_sorted_cells(file_barcode, file_bc_out, bcstart, bcend, lobc, robc)

    return bc_file_list


def clip_read(folder, read_wt):
    test = Bio.SeqIO.parse(folder + "plasmid.fasta", 'fasta')
    next(test)
    wt_gene = next(test)
    wt_gene = wt_gene.upper()

    read_wt_seq_for = read_wt
    read_wt_seq_rev = read_wt_seq_for.reverse_complement()

    if wt_gene.seq.find(read_wt_seq_for[0:18]) >= 0 or wt_gene.seq.find(read_wt_seq_for[-18:]) >= 0:
        read_wt = read_wt_seq_for
        reverse = False
    elif wt_gene.seq.find(read_wt_seq_rev[0:18]) >= 0 or wt_gene.seq.find(read_wt_seq_rev[-18:]) >= 0:
        read_wt = read_wt_seq_rev
        reverse = True
    else:
        print("could not match consensus read to plasmid sequence")

    #  Set the length of the read based on overlap with coding sequence
    if wt_gene.seq.find(read_wt) < 0:
        i = 18
        if wt_gene.seq.find(read_wt[0:18]) > 0:
            while wt_gene.seq.find(read_wt[0:i]) > 0:
                i += 1
            var_start = 1
            var_end = i - 1
        elif wt_gene.seq.find(read_wt[-18:]) > 0:
            while wt_gene.seq.find(read_wt[-i:]) > 0:
                i += 1
            var_start = -i + 2
            var_end = None
    else:
        var_start = 0
        var_end = None

    return var_start, var_end, reverse

