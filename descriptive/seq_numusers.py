"""count_users

Given a file with patterns found (results) and a sequence file, returns the
number of users matching each pattern.
"""

__version__ = '0.2'
__author__ = 'J.Renero'

import csv
import sys
import getopt
import regex as re


def checkargs(args):
    """
    Check the command line arguments and return a dictionary with all the
    input variables set.
    """
    usage = '{} -s sequence -r results_file'.format(args[0])
    args_dict = {}
    results_file = ''
    sequences_file = ''

    # We need 3 or 4 arguments
    if len(args) is not 5:
        print("Incorrect number of arguments {:d}".format(len(args)))
        print(usage)
        exit(4)

    try:
        opts, args = getopt.getopt(args[1:], "s:r:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-r"):
            results_file = arg
        elif opt in ("-s"):
            sequences_file = arg

    args_dict['results_file'] = results_file
    args_dict['sequences_file'] = sequences_file
    return args_dict


def read_results(results_file, threshold=0.0):
    """ Reads the results file (seq, #a,#b,diff,mag,f5,threshold,tg) """
    patterns = []
    with open(results_file, 'r') as input_file:
        patterns_reader = csv.reader(input_file)
        next(input_file)
        for line in patterns_reader:
            if float(line[6]) >= threshold:
                patterns.append(line[0])
    print("Read {:d} patterns.".format(len(patterns)))
    return patterns


def read_sequences(sequence_file):
    sequences = []
    with open(sequence_file, 'r') as seqs_file:
        seqs_reader = csv.reader(seqs_file)
        for line in seqs_reader:
            sequences.append(line)
    return sequences


def main(sequences_file, results_file):
    sequences = read_sequences(sequences_file)
    patterns = read_results(results_file)

    for pattern in patterns:
        print('{:40s}, '.format(pattern), end='', sep='')
        my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
        num_users = 0
        for sequence in sequences:
            if my_regex(sequence[1]):
                num_users += 1
                # print(sequence[0])
        print('{:05d}'.format(num_users))


params = checkargs(sys.argv)
main(params['sequences_file'], params['results_file'])
