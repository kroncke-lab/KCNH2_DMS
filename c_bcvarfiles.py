import pandas as pd
#  TODO: I need to keep in mind that a MAJOR feature I want is to be able to add new sequencing results directly to the already processed reads.


#  compress variant files and add counts, remove n=1 (maybe n=2, too?)
#  probably just write to csv?
def reduce_fastq(chunk, file):
    df_min = chunk.groupby(chunk.columns.tolist(), as_index=False).size()
    df_filt = df_min[df_min['size'] > 1]
    df_filt_nonwt = df_filt[df_filt['wt'] != 'wt']
    with open(file, 'a') as f:
        df_filt_nonwt.to_csv(f, header=False)

# Count unique barcode-variant pairs, discard n=1 observations (maybe also remove n=2, too?), save file 'to_feather'?
i = 0
file = folder + "data/" + str(i) + "th_fastq.tmp"
reduce_fastq(file)

df_chunks = pd.read_csv(file, names=["bc", "wt", "seq"], chunksize=2000000)

    frames = [reduce_fastq(chunk) for chunk in df_chunks]
pd.df_filt_nonwt.to_feather

