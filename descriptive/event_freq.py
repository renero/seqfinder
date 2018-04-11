"""event_freq

This module counts the frequency substrings of length 'gene_length' in
two types of files: results or patterns. 'Results' files are the output
of the 'seqfinder' process, and contains sequence strings, followed by a
number of metrics. 'Patterns' files are the input files to the 'seqfinder'
process, and contains UIDs followed by long string sequences.
The output from this module a list of the substrings found, with the
frequency of appearance in the file, and its count.
"""

__version__ = '0.2'
__author__ = 'J.Renero'


import csv
import sys
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
import getopt


def checkargs(args):
    """
    Check the command line arguments and return a dictionary with all the
    input variables set.
    """
    patterns_file = ''
    results_file = ''
    gene_length = 1
    threshold = 0.0
    usage = '{} [-r results_file | -p patterns_file] [-t threshold] [-g gene_length]'.format(args[0])
    args_dict = {}

    # We need 3 or 4 arguments
    if (3 <= len(args) <= 7) is not True:
        print("Incorrect number of arguments {:d}".format(len(args)))
        print(usage)
        exit(4)

    try:
        opts, args = getopt.getopt(args[1:], "t:r:p:g:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-r"):
            results_file = arg
        elif opt in ("-p"):
            patterns_file = arg
        elif opt in ("-t"):
            threshold = float(arg)
        elif opt in ("-g"):
            gene_length = int(arg)

    args_dict['patterns_file'] = patterns_file
    args_dict['threshold'] = threshold
    args_dict['results_file'] = results_file
    args_dict['gene_length'] = gene_length
    return args_dict


def read_results(results_file, threshold):
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


def read_patterns(patterns_file):
    """ Reads the raw sequences input file (uid, seq) """
    patterns = []
    with open(patterns_file, 'r') as input_file:
        patterns_reader = csv.reader(input_file)
        for line in patterns_reader:
            patterns.append(line[1])
    print("Read {:d} patterns.".format(len(patterns)))
    return patterns


def count_events(patterns, gene_length=1):
    counter = Counter()
    for pattern in patterns:
        if gene_length is 1:
            counter += Counter(pattern.strip())
        else:
            pat = pattern.strip()
            counter += Counter([pat[i:i+gene_length] for i in range(0, len(pat), gene_length)])
    return counter


def normalize(d, event_count, target=50.0):
    raw = sum(d.values())
    maximum_e = max(event_count, key=event_count.get)
    minimum_e = min(event_count, key=event_count.get)
    max_v = event_count[maximum_e]
    min_v = event_count[minimum_e]
    factor = target/(max_v-min_v)
    normd = {key: int(value*factor) for key, value in d.items()}
    return normd, raw


def main(results_file, patterns_file, threshold, gene_length):
    if results_file is not '':
        patterns = read_results(results_file, threshold)
    elif patterns_file is not '':
        patterns = read_patterns(patterns_file)
    event_count = count_events(patterns, gene_length)
    normalized, total = normalize(dict(event_count), event_count)
    sorted_normd = OrderedDict(sorted(normalized.items(), key=itemgetter(1), reverse=True))

    for e in sorted_normd:
        print('{}, {:.04f}, {:08d}'.format(e, (event_count[e]/total), event_count[e]))


#
# Usage: python seqstat.py [-r results_file [-t threshold]] [-p patterns_file] [-g gene_length]
#
params = checkargs(sys.argv)
main(params['results_file'], params['patterns_file'], params['threshold'], params['gene_length'])
