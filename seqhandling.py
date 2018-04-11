import logging
import collections
import csv
import statistics
import regex as re

from distsequence import DistSequence
from fuseseqs import fuse_patterns
from random import randint, random
import random
from seqlogging import setup_logging
import settings

setup_logging()
logger = logging.getLogger(__name__)


def subsequence(sequence, length_of_sequence=2, gene_length=1, max_iters=10, seed=None):
    """
    From a 'sequence' string, and an optional seed, returns a random subsequence of length 'L'.
    gene_length is the length of each event, and could be '1' for single letter seqs, or more. In case
    it is more than 1, subsequences are to be found at positions multiple of this value.
    """
    if seed is not None:
        random.seed(seed)
    if len(sequence) < length_of_sequence:
        return '_' * length_of_sequence

    valid_sequence = False
    iters = 0
    seq = ''
    while valid_sequence is False and iters < max_iters:
        random.seed()
        randpos = randint(0, len(sequence) - length_of_sequence)
        # Check if the random position of this sequence starts in the middle of a gene
        if randpos % gene_length is not 0:
            randpos -= randpos % gene_length
        seq = sequence[randpos:randpos + length_of_sequence]
        iters += 1
        if str(seq).endswith(settings.invalid_end_state) is False:
            valid_sequence = True

    return seq


def subsequences(data, num=0, seq_len=2, gene_length=1, confidence=0.05, max_search_iters=10):
    """
    Take 'num' random subsequences of length 'L' from an array of strings ('data')
    The total number of sequences are then filtered to keep only those with a prob greater
    than the confidence interval.
    'gene_length' is the length of each gene in the sequence. Default value is 1.
    'max_search_iters' is the number of times that a valid sequence is tried to be found
    by tossing random numbers, and trying to avoid that the sample ends with an invalid invalid_end_state.
    The 'invalid_end_state' is an state we want to avoid in our sequences, since it distinguishes
    sequences from one dataset perfectly from the other (like, "baja").
    """
    if num is 0:
        num = len(data) * 2
    try:
        subsequence_strings = []
        for i in range(0, num):
            valid_seq = False
            iters = 0
            subseq = ''
            while valid_seq is not True and iters < 10:
                subseq = subsequence(sequence=data[randint(0, len(data) - 1)],
                                     length_of_sequence=seq_len,
                                     gene_length=gene_length,
                                     max_iters=max_search_iters)
                if str(subseq).endswith(settings.invalid_end_state) is False:
                    valid_seq = True
                iters += 1
            subsequence_strings.append(subseq)
    except TypeError:  # Some parts of data may contain null or nan and simply fail.
        logger.error("TypeError on len(data)=", len(data), ", num=", num, ", and L=", seq_len)
        subsequence_strings = []
    distribution = collections.Counter(subsequence_strings)
    relevant_subseqs = [key for key in distribution.keys() if (distribution[key] / num) > confidence]
    return list(set(relevant_subseqs))


def avg_count_pattern_in_data(clusters, data):
    """
    Count the number of occurrences of the prototypes in clusters in the data.
    """
    prototypes = []
    for cluster in clusters:
        prototypes.append(cluster.prototype)
    pattern_counts = pattern_count_seqs(data, prototypes)
    avg_count = statistics.mean([pattern_count.count for pattern_count in pattern_counts])
    return avg_count


def howmany(my_regex, string):
    """
    Returns how many occurrences of the pattern can be found in the string.
    """
    try:
        matches = my_regex(string)
    except TypeError:  # Some strings are null for no known reason
        matches = []
        logger.error("TypeError on pattern:", my_regex, "and string:", string)
    return len(matches)


def pattern_count_kprototype(data, k_prototype):
    """
    Return the number of overlapped occurrences of pattern in data
    - Data is an array of strings, where the patterns is searched for.
    - kPrototype is KPrototype object.
    """
    pattern_count = []
    for ptype in k_prototype.prototype:
        pattern = ptype.prototype
        pat_count = 0
        my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
        for string in data:
            pat_count += howmany(my_regex, string)
        pattern_count.append((pattern, pat_count))
    uniq_tuples = list(set(pattern_count))
    return sorted(uniq_tuples, key=lambda x: x[1], reverse=True)


def pattern_count_prototype(data, prototypes):
    """
    Return the number of overlapped occurrences of pattern in data
    - Data is an array of strings, where the patterns is searched for.
    - Protoypes is an array of Prototype objects.
    """
    logger.debug("Computing prototypes count in data...")
    pattern_count = []
    for ptype in prototypes:
        pattern = ptype.prototype
        pat_count = 0
        my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
        for string in data:
            pat_count += howmany(my_regex, string)
        pattern_count.append((pattern, pat_count))
    uniq_tuples = list(set(pattern_count))
    return sorted(uniq_tuples, key=lambda x: x[1], reverse=True)


def pattern_count_seqs(data, patterns_seqs):
    """
    Return the number of overlapped occurrences of pattern in data
    - Data is an array of strings, where the patterns is searched for.
    - Pattern_seqs is an array of strings representing the sequence prototypes.
    """
    from collections import namedtuple
    PatternCount = namedtuple('PatternCount', ['pattern', 'count'])
    pattern_count = []
    for pattern in patterns_seqs:
        pat_count = 0
        my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
        for string in data:
            pat_count += howmany(my_regex, string)
        pattern_count.append(PatternCount(pattern, pat_count))
    uniq_tuples = list(set(pattern_count))
    return sorted(uniq_tuples, key=lambda x: x.count, reverse=True)


def distinctive_sequences(candidate_prototypes, data_bajas, data_altas, filter_distinctive=False,
                          minimum_gscore=0.5):
    """
    Take the prototypes from 'bajas' and count them in 'altas'. From the entire list of prototypes
    filter only those in altas in less than 'max_ocurrences'.
    G_Score is the ratio between count in bajas and count in bajas+altas.
    :param filter_distinctive:
    :return: Distinctive Sequence composed of patterns,countInA, countInB
    """

    # Generate also the fused patterns
    logger.info("Original patterns list size: %d", len(candidate_prototypes))
    candidate_prototypes.extend(fuse_patterns(candidate_prototypes))
    logger.info("Extended with fused patterns size: %d", len(candidate_prototypes))
    candidate_prototype_count_altas = pattern_count_seqs(data_altas, candidate_prototypes)

    # Determine the max nr of occurrences of prototypes in 'altas'. If a pattern is present more than
    # this max_counts, then it is discarded.
    sum_counts = sum([prototype.count for prototype in candidate_prototype_count_altas])
    logger.debug("Candidate prototypes count in B-set (altas): %d", sum_counts)

    # The the maximum number of times that I allow each sequence to appear in the other dataset tp 5% of the total
    # number of occurrences. If that number is not possible, set it to 1% of the total length. And if that is
    # not possible either, set it to the length of altas
    if filter_distinctive is True:
        max_ocurrences = int(sum_counts * 0.05)
        if max_ocurrences < 10:
            max_ocurrences = int(len(data_altas) * 0.01)
            if max_ocurrences is 0:
                max_ocurrences = len(data_altas)
        logger.debug("Setting threshold for occurrences in B-Set to: %d", max_ocurrences)
        rare_candidate_patterns = list(filter(lambda x: x.count < max_ocurrences, candidate_prototype_count_altas))
    else:
        rare_candidate_patterns = list(candidate_prototype_count_altas)

    rare_altas_patterns_strings = [rare_cand_pat.pattern for rare_cand_pat in rare_candidate_patterns]
    rare_candidate_patterns_count_bajas = pattern_count_seqs(data_bajas, rare_altas_patterns_strings)

    # Join the two lists.
    dict1 = dict(rare_candidate_patterns)
    dict2 = dict(rare_candidate_patterns_count_bajas)

    lst3 = [(k, dict1[k], dict2[k]) for k in sorted(dict1)]
    dist_seqs = list(filter(lambda x: x.diff > 0,
                            [DistSequence(lst3_elem[0], lst3_elem[1], lst3_elem[2]) for lst3_elem in lst3]))

    # return only those sequences with g_score greater or equal than the minimum specified.
    filtered_list = list(filter(lambda ds: ds.gscore >= minimum_gscore, dist_seqs))
    logger.debug("Filtered list size: %d", len(filtered_list))
    return filtered_list


def show_dss(dss):
    dss = list(sorted(dss, key=lambda ds: ds.gscore, reverse=True))
    tp = 0
    tn = 0
    print("")
    print("  Pattern            | # Altas | # Bajas |   Diff |Magnitude | Goodness |  TPR  |  G·TPR  ")
    print("  -------------------|---------|---------|--------|----------|----------|-------|---------")
    for x in dss:
        print("{:>20} | {:7d} | {:7d} | {:6d} | {:8.1f} | {:8.1f} | {:5.2f} | {:7.2f}"
              .format(x.pattern, x.a_count, x.b_count, x.diff, x.ratio, x.goodness, x.tpr, x.gscore))
        tp += x.b_count
        tn += x.a_count
    print("")
    print("  T.Bajas = {:d}, T.Altas = {:d}".format(tp, tn))
    if tp + tn is not 0:
        print("  TPR......= {:5.2f}".format(tp / (tp + tn)))
    else:
        print("  TPR......= {:5.2f}".format(0))


def save_dss(dss, outpath):
    dss = list(sorted(dss, key=lambda ds: ds.gscore, reverse=True))
    writer = csv.writer(open(outpath, 'w'))
    writer.writerow(["pattern", "altas", "bajas", "diff", "ratio", "goodness", "tpr", "gtpr"])
    for x in dss:
        writer.writerow([x.pattern, x.a_count, x.b_count, x.diff, x.ratio, x.goodness, x.tpr, x.gscore])
