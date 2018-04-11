import math


class DistSequence(object):
    def __init__(self, pattern, count_in_a, count_in_b):
        self.pattern = pattern
        self.a_count = count_in_a
        self.b_count = count_in_b
        self.diff = self.b_count - self.a_count
        self.ratio = (self.b_count / (self.a_count + 0.1))
        self.goodness = (self.diff * (self.ratio / 100)) / math.log2(self.a_count + 2)
        if (self.a_count + self.b_count) is not 0:
            self.tpr = self.b_count / (self.a_count + self.b_count)
        else:
            self.tpr = 0
        self.gscore = self.goodness * self.tpr
