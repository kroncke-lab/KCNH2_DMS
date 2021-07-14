from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re
import models.consensus
import matplotlib.pyplot as plt
import pandas as pd

def check_align(sequence, data, index):
    unmatched = 0
    for line in data:
        if not re.search(sequence[index:index + 18], line):
            unmatched = unmatched + 1
    return unmatched


file_in = files_barcode[1]
consensus_seq = models.consensus.findconsensus(file_in)
data = []
i = 0
unmatched = 0
with open(file_in) as in_handle:
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle):
        print(seq)
        data.append(seq)
        count += 1
        if count == 10000:
            break

appended_data = []
prev_match = 0
diff = 0
cumulative = 0
barcode = False
for ind in range(1, 60, 1):
    match = check_align(bc, data, ind)
    if prev_match > 0:
        diff = match - prev_match
    prev_match = match
    row = [ind, match, diff]
    appended_data.append(row)
    cumulative = cumulative + match
    if cumulative > 180000:
        barcode = True

df = pd.DataFrame(appended_data, columns=["index", "matches", "diff"])
df.plot.line()
lower = df[['diff']].idxmax()
upper = df[['diff']].idxmin()
lobc = bc[upper.values[0]-18-5:upper.values[0]-18-1]
robc = bc[upper.values[0]+1:upper.values[0]+5]

for d in data:
    ds = d[upper.values[0]-18-7:upper.values[0]+7]
    print(ds)
    m = re.split(lobc, ds[0:7])
    print(len(m[0]))
    print(m)
    ds = ds[len(m[0])+len(lobc):len(ds)]
    print(ds)
    try:
        if m[1]:
            n = re.split(robc, ds[17:len(ds)])
            print(len(n[0]))
            print(ds[0:len(ds)-len(n[0])])
    except IndexError:
        continue


#  I can assume the maximum shift will be +/- 1 for nearly all reads
#  so a reasonable estimate will be upper - 18
#  then I can use (upper - 18) - 5 and upper + 5 as the indices to pull
#  from each sequence, then grab upper - 18 - 5 to upper - 18 - 1 (one side of barcode)
#  and upper + 1 to upper + 4

tmp = []
data = []
appended_data = []
prev_match = 0
diff = 0
cumulative = 0
barcode = False
with open(files[0]) as in_handle:
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle):
        tmp.append(Seq(seq, IUPAC.ambiguous_dna))
        data.append(seq)
        count += 1
        if count == 40000:
            wt = motifs.create(tmp, IUPAC.ambiguous_dna).consensus
            break
print(wt._data)

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
