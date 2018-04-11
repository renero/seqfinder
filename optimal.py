import statistics

from patterns import find_patterns
from kprototype import KPrototype
from seqhandling import avg_count_pattern_in_data


def find_elbow(curve):
    import numpy as np

    """
    Find the elbow of curve shaped like the inverse of a log.
    Taken from: https://stackoverflow.com/a/37121355/892904
    """
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    np.array([range(nPoints), curve])
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec ** 2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(vecFromFirst * np.tile(lineVecNorm, (nPoints, 1)), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint


def best_prototype(clustering, Kprototypes, compute_f3=False):
    """
    Returns the best prototype form the many different ones generated with the differente clusterings
    """
    # Best clustering according to cluster intra distance behavior
    converg = [clustering[i].eq for i in range(0, len(clustering) - 2)]
    if converg is not None:
        f0 = converg.index(min(converg))
    else:
        f0 = 0

    # Best Clustering Prototypes according to the matching rate.
    norm_rates = [kprototype.mean_rate / len(Kprototypes) for kprototype in Kprototypes]
    f1 = norm_rates.index(max(norm_rates))

    f2 = find_elbow(norm_rates)
    if f0 >= len(clustering) / 2:
        f0 = f2
    # print("Elbow(%d)=%01.2f, F0(%d)=%01.2f, F1(%d)=%01.2f" %
    #    (elbow,Kprototypes[elbow].fscore,f0,Kprototypes[f0].fscore,f1,Kprototypes[f1].fscore))

    #  Best prototype as the one with minimum number of occurrences of its prototypes in another dataset.
    if compute_f3:
        avgs = []
        for kprototype in Kprototypes:
            avgs.append(kprototype.avg_count)
        f3 = avgs.index(min(avgs))
    else:
        f3 = 0

    return f0, f1, f2, f3


def obtain_kprototypes(clustering, regex_length=None, ctrl_data=None, gene_length=1):
    """
    Obtain the prototypes for each cluster on each of the different clustering results obtained.
    Returns a 3-tuple with the prototypes, the mean score and the mean rate for each clustering set.
    :param length:
    """
    k_prototypes = []
    for cluster_group in clustering:
        clusters = [find_patterns(cluster_group.clusters[i], regex_length, gene_length=gene_length)
                    for i in cluster_group.clusters.keys()]
        scores = [clusters[i].score for i in range(0, len(clusters))]
        rates = [clusters[i].rate for i in range(0, len(clusters))]
        if ctrl_data is not None:
            avg_count = avg_count_pattern_in_data(clusters, ctrl_data)
        else:
            avg_count = 0
        k_prototypes.append(KPrototype(clusters, statistics.mean(scores), statistics.mean(rates), avg_count))

    return k_prototypes
