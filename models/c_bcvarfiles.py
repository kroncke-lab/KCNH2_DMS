import pandas as pd
from functools import reduce

# TODO: I need to keep in mind that a MAJOR feature I want is to be able to add new sequencing results directly
#  to the already processed reads.


def get_counts(chunk):
    uniq_bc_pair = chunk[["bc", "seq"]]
    return uniq_bc_pair.value_counts()


def get_counts_sort(chunk):
    uniq_bc_pair = chunk[["bc"]]
    return uniq_bc_pair.value_counts()


def add(previous_result, new_result):
    print("Added a new chunk! (500,000 records)")
    return previous_result.add(new_result, fill_value=0)


def find_and_reduce_bcvars(file):
    chunks = pd.read_csv(file, names=["bc", "seq"], chunksize=500000)
    processed_chunks = map(get_counts, chunks)
    #  this next part takes ~ 15 min.
    result = reduce(add, processed_chunks)
    # la siquiente linea dejan de los barcode-variantes llaves visto solo uno o dos veces.
    result = result[result > 4]
    return result


def find_and_reduce_bc_sort(file):
    chunks = pd.read_csv(file, names=["bc"], chunksize=500000)
    processed_chunks = map(get_counts_sort, chunks)
    result = reduce(add, processed_chunks)
    return result

