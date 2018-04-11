"""count_users

Search a given pattern in a sequence file, and return the UID of
all users matching it.
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
    usage = '{} [-c | -l] -p pattern -s sequence_file'.format(args[0])
    args_dict = {}
    sequences_file = ''
    pattern = ''
    count_users = False
    list_users = False

    # We need 3 or 4 arguments
    if len(args) is not 6:
        print("Incorrect number of arguments {:d}".format(len(args)))
        print(usage)
        exit(4)

    try:
        opts, args = getopt.getopt(args[1:], "cls:p:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-s"):
            sequences_file = arg
        elif opt in ('-c'):
            count_users = True
        elif opt in ('-l'):
            list_users = True
        elif opt in ("-p"):
            pattern = arg

    args_dict['sequences_file'] = sequences_file
    args_dict['pattern'] = pattern
    args_dict['count_users'] = count_users
    args_dict['list_users'] = list_users
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


def main(pattern, sequences_file, count_em, list_em):
    sequences = read_sequences(sequences_file)
    print('{:s}'.format(pattern), end='', sep='')
    if list_em is True:
        print('')
    else:
        print(',', end='')
    my_regex = re.compile(r'((' + pattern.replace('?', '\w') + '))').findall
    num_users = 0
    for sequence in sequences:
        if my_regex(sequence[1]):
            num_users += 1
            if list_em is True:
                print(sequence[0])
    if count_em is True:
        print('{:d}'.format(num_users), sep='')


params = checkargs(sys.argv)
main(params['pattern'], params['sequences_file'], params['count_users'],
     params['list_users'])
