import logging
from seqlogging import setup_logging

from cluster import clusterize
from optimal import obtain_kprototypes, best_prototype
from pseq import PSeq
from seqhandling import subsequences, distinctive_sequences

setup_logging()
logger = logging.getLogger(__name__)


def seqfind(data, ctrl_data=None, num_samples=100, seq_length=3, regex_length=None, seq_min_count_percentage=0.05,
            compute_f3=False, exhaustive_search_mode=False, gene_length=1):
    """
    Given a dataset ('data'), returns the different clusterings for the relevant subsequences found
    on it of length 'seq_length'. Once the different clusterings have been produced, the best prototypes
    for each clustering are computed according to three different criteria: intra-cluster distance (clusters
    with homogeneus strings are favoured), max score in matching (prototypes matching the larger number of
    individuals within each cluster), and 'elbow' which is the point where the second derivative of the curve
    representing the different matching scores reaches its maximum point. This last method is the same used
    to compute the optimal value in ROC curves.

    Arguments:
    - data: The dataset as an array of strings
    - num_samples: The number of samples taken from the dataset to start working with.
    - seq_length: The length of each sequence extracted from the strings in data.
    - regex_length: The maximum nr of wildcards to be used when searching for patterns.
    - draw: Boolean indicating if summary drawings must be produced or not.

    Returns:
    - The clustering
    - The prototypes for each clustering
    - The values for optimal prototypes found according the three different metrics.
    :param gene_length:
    """
    subseqs = subsequences(data, num=num_samples, seq_len=seq_length, confidence=seq_min_count_percentage,
                           gene_length=gene_length)
    subseqs_ctrl = subsequences(ctrl_data, num=num_samples, seq_len=seq_length, confidence=seq_min_count_percentage,
                                gene_length=gene_length)
    clustering = [clusterize(subseqs, num_clusters=k, gene_length=gene_length) for k in range(2, len(subseqs) - 1)]
    kprototypes = obtain_kprototypes(clustering, regex_length, subseqs_ctrl, gene_length)

    if len(clustering) == 0:
        logger.error("ERROR: No clustering possible. Try decreasing min count percentage for sequence patterns.")
        return None

    if exhaustive_search_mode is False:
        # Best Clustering Prototypes according to the different metrics
        f0, f1, f2, f3 = best_prototype(clustering, kprototypes, compute_f3)
    else:
        f0 = f1 = f2 = f3 = 0
    return PSeq(clustering, kprototypes, f0, f1, f2, f3)


def take_only_newones(candidates, final_list):
    if final_list is []:
        return candidates
    s1 = set(candidates)
    s2 = set([element.pattern for element in final_list])
    return list(s1 - s2)


def best_distinctive(ds_list, data_bajas, data_altas, num_samples=100, seq_length=3, regex_length=None, ntimes=10,
                     seq_min_count_percentage=0.05, min_g_score=0.5, filter_distinctive=True,
                     exhaustive_search=False, gene_length=1):
    for n in range(ntimes):
        logger.info("  Iteration {:02d}/{:02d}".format(n + 1, ntimes))
        seqs = seqfind(data_bajas, data_altas,
                       num_samples=num_samples,
                       seq_length=seq_length,
                       regex_length=regex_length,
                       seq_min_count_percentage=seq_min_count_percentage,
                       compute_f3=True,
                       exhaustive_search_mode=exhaustive_search,
                       gene_length=gene_length)

        # Search for good patterns in the candidates selected by the 4 metrics
        if exhaustive_search:
            if seqs is not None:
                search_range = range(len(seqs.kprototypes))
            else:
                logger.warning("No prototypes selected.")
                search_range = []
        else:
            search_range = (seqs.f0, seqs.f1, seqs.f2, seqs.f3)
        candidate_prototypes = []
        for cluster_index in search_range:
            candidate_prototypes = take_only_newones(seqs.kprototypes[cluster_index].prototypes(), ds_list)

        # Search only those really distinctive. Control relevance with GScore and the filtering flag.
        ds = distinctive_sequences(candidate_prototypes, data_bajas, data_altas,
                                   filter_distinctive=filter_distinctive, minimum_gscore=min_g_score)
        logger.debug('Added {:d} signatures to global list.'.format(len(ds)))
        ds_list.extend(ds)
    if len(ds_list) is not 0:
        logger.debug("Final list of {:d} signatures".format(len(ds_list)))
