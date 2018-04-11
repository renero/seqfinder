import logging
from seqlogging import setup_logging
import sys

from seqhandling import pattern_count_prototype


setup_logging()
logger = logging.getLogger(__name__)

class PSeq(object):
    def __init__(self, clustering=None, kprototypes=None, f0=0.0, f1=0.0, f2=0.0, f3=0.0):
        self.clustering = clustering
        self.kprototypes = kprototypes
        self.f0 = f0
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3

    def best_score(self):
        if self.kprototypes is None:
            return 0.0
        scores = [self.kprototypes[self.f0].fscore,
                  self.kprototypes[self.f1].fscore,
                  self.kprototypes[self.f2].fscore]
        return max(scores)

    def best_idx(self):
        """Returns the index of the prototype with larger with larger score (F0, F1 or F2)"""
        if self.kprototypes is None:
            return
        scores = [self.kprototypes[self.f0].fscore,
                  self.kprototypes[self.f1].fscore,
                  self.kprototypes[self.f2].fscore]
        best = max(scores)
        return [self.f0, self.f1, self.f2, self.f3][scores.index(best)]

    def show_f0(self):
        self.kprototypes[self.f0].show()

    def show_f1(self):
        self.kprototypes[self.f1].show()

    def show_f2(self):
        self.kprototypes[self.f2].show()

    def show_f3(self):
        self.kprototypes[self.f3].show()

    def show_proto(self, idx):
        self.kprototypes[idx].show()

    def show_best(self):
        if self.kprototypes is None:
            return
        idx = self.best_idx()
        self.kprototypes[idx].show()

    def best_prototype_str(self, data):
        """Returns the string of the best prototype in terms of occurrences in the data."""
        if self.kprototypes is None:
            return
        idx = self.best_idx()
        pc = pattern_count_prototype(data, self.kprototypes[idx].prototype)
        return pc[0][0]

    def best_prototypes_str(self, data, count, threshold):
        """
        Returns an array (count elements) with the prototypes beyond a threshold
        (# of occurrences in the data)
        """
        if self.kprototypes is None:
            logger.error("ERROR: Empty KPrototypes")
            return
        idx = self.best_idx()
        pc = pattern_count_prototype(data, self.kprototypes[idx].prototype)
        pc = list(filter(lambda x: x[1] > threshold, pc))
        return [pc[i][0] for i in range(min([count, len(pc)]))]

    def print_kprototypes(self, idx, data, threshold):
        if self.kprototypes is None:
            logger.info("ERROR: Empty KPrototypes")
            return
        pc = pattern_count_prototype(data, self.kprototypes[idx].clusters)
        pc = list(filter(lambda x: x[1] > threshold, pc))
        sys.stdout.flush()
        logger.info(end='')
        logger.info('Sequence        Count')
        logger.info('---------------------')
        for pat_count in pc:
            logger.info('{:15}{:>5d}'.format(pat_count[0], pat_count[1]))
