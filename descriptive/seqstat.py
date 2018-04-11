#!/usr/bin/env python
"""seqstat

This module produces basic statistics for the sequences found in input file.
"""

__version__ = '0.1'
__author__ = 'J.Renero'


import os,sys
import csv
import numpy as np


def read_patterns(patterns_file):
    """ Reads the raw sequences input file (uid, seq) """
    patterns = []
    with open(patterns_file, 'r') as input_file:
        patterns_reader = csv.reader(input_file)
        for line in patterns_reader:
            patterns.append(line[1])
    print("count,{:d}".format(len(patterns)))
    return patterns


def basic_stats(input_file):
    patterns = read_patterns(input_file)
    lengths = np.array([])
    for pattern in patterns:
        lengths = np.append(lengths, len(pattern))
    print('mean,{:.04}'.format(np.mean(lengths)))
    print('median,{:.04}'.format(np.median(lengths)))
    hist, bin_edges = np.histogram(lengths, bins=30, range=(0.0,300.0))
    print('histogram,', ','.join([str(h) for h in hist]), sep='')
    print('hst_edges,', ','.join([str(int(bh)) for bh in bin_edges]), sep='')


def main(argv):
    if len(sys.argv) is not 2:
        print("Usage: {} <patterns_file>".format(sys.argv[0]))
        exit(1)
    else:
        input_file = sys.argv[1]

    basic_stats(input_file)


main(sys.argv)
