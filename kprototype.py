import logging
from seqlogging import setup_logging

setup_logging()
logger = logging.getLogger(__name__)


class KPrototype(object):
    def __init__(self, cluster_list, mean_score, mean_rate, avg_count):
        self.clusters = cluster_list
        self.mean_rate = mean_rate
        self.fscore = mean_score * mean_rate
        self.norm_rate = mean_rate / len(cluster_list)
        self.avg_count = avg_count  # Avg nr. of occurrences of patterns in control data.

    def size(self):
        return len(self.clusters)

    def show(self):
        logger.info("K = %d. Avg.Rate: %02.2f, F-Score=%1.02f, Avg Count=%d" %
              (len(self.clusters), self.mean_rate, self.fscore, self.avg_count))
        for p in (sorted(self.clusters, key=lambda x: x.rate, reverse=True)):
            p.show()

    def best(self):
        return sorted(self.clusters, key=lambda x: x.rate, reverse=True)[0].prototype

    def prototypes(self):
        return [cluster.prototype for cluster in self.clusters]
