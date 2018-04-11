"""count_users

This module receives as inputs a list of patterns (results) and two files with uids and patterns.
All patterns in results are searched in both files (positive and negative) and the number of lines
where the pattern is found is returned and summed up.
A threshold in the results file can also be used to control which patterns to search in both
files. Therefore, for each threshold a true-positive-rate and true-negative-rate are returned
corresponding to the number of times patterns in results file, above a given threshold are found
on each set.
"""

__version__ = '0.2'
__author__ = 'J.Renero'

import csv
import numpy as np
import sys
import getopt
import regex as re


def checkargs(args):
    """
    Check the command line arguments and return a dictionary with all the
    input variables set.
    """
    bajas_file = ''
    altas_file = ''
    results_file = ''
    gene_length = 1
    threshold_ini = 0.0
    threshold_end = 0.0
    usage = '{} -r results_file -b bajas_file -a altas_file -i threshold_ini -e threshold_end'.format(args[0])
    args_dict = {}

    # We need 3 or 4 arguments
    if (6 <= len(args) <= 11) is not True:
        print("Incorrect number of arguments {:d}".format(len(args)))
        print(usage)
        exit(4)

    try:
        opts, args = getopt.getopt(args[1:], "r:b:a:i:e:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-r"):
            results_file = arg
        elif opt in ("-b"):
            bajas_file = arg
        elif opt in ("-a"):
            altas_file = arg
        elif opt in ("-i"):
            threshold_ini = float(arg)
        elif opt in ("-e"):
            threshold_end = float(arg)

    args_dict['bajas_file'] = bajas_file
    args_dict['altas_file'] = altas_file
    args_dict['threshold_ini'] = threshold_ini
    args_dict['threshold_end'] = threshold_end
    args_dict['results_file'] = results_file
    return args_dict


def read_results(results_file, threshold):
    patterns = []
    with open(results_file, 'r') as patterns_file:
        patterns_reader = csv.reader(patterns_file)
        next(patterns_file)
        for line in patterns_reader:
            if float(line[6]) >= threshold:
                patterns.append(line[0])
    return patterns


def read_patterns(sequence_file):
    sequences = []
    with open(sequence_file, 'r') as seqs_file:
        seqs_reader = csv.reader(seqs_file)
        for line in seqs_reader:
            sequences.append(line)
    return sequences


def find_users_matching_patterns(sequence_file, patterns):
    sequence = read_patterns(sequence_file)
    # This named tuple holds the different matches found in the dataset
    uids_matching = []
    for sequence in sequences:
        for pattern in patterns:
            if pattern in sequence[1]:
                uids_matching.append(sequence[0])
                break
    return list(set(uids_matching))


def find_users_per_sequence(results_file, positive_file, negative_file):
    patterns = read_results(results_file, 0.0)
    positive_users = read_patterns(positive_file)
    negative_users = read_patterns(negative_file)
    for pattern in patterns:
        positive_matches = 0
        negative_matches = 0
        my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
        for uid in positive_users:
            if my_regex(uid[1]):
                positive_matches += 1
        for uid in negative_users:
            if my_regex(uid[1]):
                negative_matches += 1
        print('{:s},{:d},{:d}'.format(pattern, positive_matches, negative_matches))


def main(results_file, positive_sequence_file, negative_sequence_file, threshold_ini, threshold_end):
    find_users_per_sequence(results_file, positive_sequence_file, negative_sequence_file)
    print('{:<10s};{:<10s};{:<10s};{:<10s}'.format('threshold', 'patterns', 'positives', 'negatives'))
    for threshold in np.arange(threshold_ini, threshold_end, 0.01):
        patterns = read_results(results_file, threshold)
        positives = find_users_matching_patterns(positive_sequence_file, patterns)
        negatives = find_users_matching_patterns(negative_sequence_file, patterns)
        print('{:<10.2f};{:<10d};{:<10d};{:<10d}'.format(threshold, len(patterns), len(positives), len(negatives)))


params = checkargs(sys.argv)
main(params['results_file'], params['bajas_file'], params['altas_file'], 
    params['threshold_ini'], params['threshold_end'])
