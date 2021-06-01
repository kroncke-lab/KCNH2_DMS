import pandas as pd
from functools import reduce

# TODO: I need to keep in mind that a MAJOR feature I want is to be able to add new sequencing results directly
#  to the already processed reads.


def get_counts(chunk):
    uniq_bc_pair = chunk[["bc", "wt", "seq"]]
    return uniq_bc_pair.value_counts()


def add(previous_result, new_result):
    return previous_result.add(new_result, fill_value=0)


def find_and_reduce_bcvars(file):
    chunks = pd.read_csv(file, names=["bc", "wt", "seq"], chunksize=500000)
    processed_chunks = map(get_counts, chunks)
    #  this next part takes ~ 15 min.
    result = reduce(add, processed_chunks)
    # la siquiente linea deja de los barcode-variadades llaves visto solo uno o dos veces.
    result = result[result > 5]
    return result


