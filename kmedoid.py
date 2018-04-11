import logging
from seqlogging import setup_logging
import statistics


setup_logging()
logger = logging.getLogger(__name__)


class KMedoid(object):
    def __init__(self, clusters, intra_distance):
        self.clusters = clusters
        self.intra_distance = intra_distance
        self.mean = statistics.mean(intra_distance)
        self.median = statistics.median(intra_distance)
        self.stdev = statistics.stdev(intra_distance)
        self.variance = statistics.variance(intra_distance)
        # The point between the |avg(mean,median) - stdev| representing
        # the best possible combination of the three.
        self.eq = abs(statistics.mean([self.mean, self.median]) - self.stdev)

    def summary(self):
        logger.info("K = %d. Mean: %02.2f, Median: %02.2f, Stdev: %02.2f" %
                    (len(self.clusters), self.mean, self.median, self.stdev))
        counter = 0
        for label, grp in self.clusters.items():
            logger.info("Cluster %02d: %s" % (counter, ','.join(grp)))
            logger.info(" intr_dist: %02.2f" % self.intra_distance[counter])
            counter += 1
