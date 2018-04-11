import numpy as np


def init():
    global seed, min_g_score, filter_distinctive
    global hamming_vectorized, levehnstein_vectorized, distance_function
    global NOI, LOS, HMW, MPT, MCP, NOS
    global exhaustive_search, invalid_end_state, gene_length
    global bajas, altas, min_seq_len, max_seq_len, delimiter

    bajas = '~/Code/seqfinder/data3/seqs_bajas_ncnf_norep_v3.1.csv'
    altas = '~/Code/seqfinder/data3/seqs_altas_ncnf_norep_v3.1.csv'
    delimiter = ','

    seed = 4
    gene_length = 2
    min_seq_len = 6
    max_seq_len = 8

    hamming_vectorized = np.vectorize(hamming_distance)
    levehnstein_vectorized = np.vectorize(levenhstein_distance)
    distance_function = hamming_vectorized

    exhaustive_search = False
    filter_distinctive = False
    invalid_end_state = 'BB'                     # set to empty if no end state to avoid.

    # Number-of-Iterations
    NOI = 10
    # How many wildcards accepted?
    HMW = 1
    # Threshold for the nr. of occurrences in data of the prototypes found.
    MPT = 1
    # Minimum count percentage of presence of each sequence in the data to be representative.
    MCP = 0.00000
    # Number of samples to take each time
    NOS = 100
    # Minimum ratio between sequence presence COUNT in bajas and altas.
    # This is the so called G-Score and is #InBajas/#InBajas+InAltas.
    # If 1.0 means that the sequence must not be present in Altas.
    # If 0.0 means that the sequence is more present in Altas than in Bajas.
    min_g_score = 0.9
