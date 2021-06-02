import pandas as pd
#  TODO: I need to keep in mind that a MAJOR feature I want is to be able to add new sequencing results directly to the already processed reads.
#  TODO: I need to compare each line with the "wt" discovered in the first part of the script.


#  compress variant files and add counts, remove n=1 (maybe n=2, too?)
#  probably just write to csv?
def reduce_fastq(chunk, file):
    df_min = chunk.groupby(chunk.columns.tolist(), as_index=False).size()
    df_filt = df_min[df_min['size'] > 1]
    df_filt_nonwt = df_filt[df_filt['wt'] != 'wt']
    with open(file, 'a') as f:
        df_filt_nonwt.to_csv(f, header=False)

# Count unique barcode-variant pairs, discard n=1 observations (maybe also remove n=2, too?), save file 'to_feather'?
folder = "E:\\Tile1-Library2-LV310\\barcode-key\\20201006-5406-LVkeyRepeat\\"  # input file

i = 0
file = folder + "data/" + str(i) + "th_fastq.tmp"
reduce_fastq(file)

df_chunks = pd.read_csv(file, names=["bc", "wt", "seq"], chunksize=2000000)

    frames = [reduce_fastq(chunk) for chunk in df_chunks]
pd.df_filt_nonwt.to_feather


from functools import reduce

def get_counts(chunk):
    voters_street = chunk[
        "Residential Address Street Name "]
    return voters_street.value_counts()

def add(previous_result, new_result):
    return previous_result.add(new_result, fill_value=0)

# MapReduce structure:
chunks = pandas.read_csv("voters.csv", chunksize=1000)
processed_chunks = map(get_counts, chunks)
result = reduce(add, processed_chunks)

result.sort_values(ascending=False, inplace=True)
print(result)
