import distance
import logging
import numpy as np
from itertools import zip_longest
from kmedoid import KMedoid


import settings
from seqlogging import setup_logging

setup_logging()
logger = logging.getLogger(__name__)


def kMedoids(D, k, tmax=100):
    # determine dimensions of distance matrix D
    m, n = D.shape

    if k > n:
        raise Exception('too many medoids')
    # randomly initialize an array of k medoid indices
    M = np.arange(n)
    np.random.shuffle(M)
    M = np.sort(M[:k])

    # create a copy of the array of medoid indices
    Mnew = np.copy(M)

    # initialize a dictionary to represent clusters
    C = {}
    for t in range(tmax):
        # determine clusters, i. e. arrays of data indices
        J = np.argmin(D[:, M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J == kappa)[0]
        # update cluster medoids
        for kappa in range(k):
            J = np.mean(D[np.ix_(C[kappa], C[kappa])], axis=1)
            j = np.argmin(J)
            Mnew[kappa] = C[kappa][j]
        np.sort(Mnew)
        # check for convergence
        if np.array_equal(M, Mnew):
            break
        M = np.copy(Mnew)
    else:
        # final update of cluster memberships
        J = np.argmin(D[:, M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J == kappa)[0]

    # return results
    return M, C


def seq_strip(sequence, gene_length=1, fillvalue='_'):
    """Collect data into fixed-length chunks or blocks"""
    # seq_strip('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    if gene_length is 1:
        return sequence
    args = [iter(sequence)] * gene_length
    zl_object = zip_longest(*args, fillvalue=fillvalue)
    return [''.join(x) for x in list(zl_object)]


def hamming_distance(s1, s2, gene_length):
    return distance.hamming(seq_strip(s1, gene_length), seq_strip(s2, gene_length), normalized=True)


def levenhstein_distance(s1, s2, gene_length):
    return distance.levenshtein(seq_strip(s1, gene_length), seq_strip(s2, gene_length), normalized=True)


def distance_matrix(subseqs, gene_length=1, triangle=True):
    dist_matrix = settings.distance_function(subseqs[:, np.newaxis], subseqs, gene_length)
    if triangle is True:
        indices = np.tril_indices_from(dist_matrix, -1)
        return dist_matrix[indices]
    else:
        return dist_matrix


def avg_intracluster_distance(clusters):
    avg = []
    for label, grp in clusters.items():
        cluster_elements = grp
        d = distance_matrix(np.array(cluster_elements), triangle=True)
        if len(d) > 0:
            avg.append(sum(d) / len(d))
        else:
            avg.append(sum(d))
    return avg


def clusterize(subseqs, num_clusters, gene_length=1):
    """
    Run kMedoids over the list of unique subsequences using 'k' clusters.
    Distance function can be modified using 'method' argument as 'L' for Levenshtein or 'H' for Hamming.
    """
    # Set 'triangle' to True if using BioPython KMedoids library.
    dist = distance_matrix(np.array(subseqs), gene_length=gene_length, triangle=False)
    logger.debug("Calling KMedoids with list length=%d and num_cluster=%d", len(dist.tolist()), num_clusters)
    try:
        #labels, error, nfound = PC.kmedoids(dist.tolist(), nclusters=num_clusters)
        _, labels = kMedoids(dist, num_clusters)
    except RuntimeError as err:
        logger.warning("Runtime error when trying to launch KMedoid: %s", err)
        logger.warning("dist.tolist() = %s", str(dist.tolist()))
        logger.warning("nclusters = %d", num_clusters)
    cluster = dict()
    for label in labels:
        [cluster.setdefault(label, []).append(subseqs[idx]) for idx in labels[label]]
    # Uncomment below if using BioPython kMedoids.
    # for word, label in zip(subseqs, labels):
    #     cluster.setdefault(label, []).append(word)
    intra_distance = avg_intracluster_distance(cluster)
    return KMedoid(cluster, intra_distance)
