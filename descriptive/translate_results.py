import os, sys
import csv

def read_dict(dict_file):
    gene_dict = dict()
    with open(dict_file, 'r') as dictionary_file:
        dictionary_reader = csv.reader(dictionary_file)
        for line in dictionary_reader:
            gene_dict.setdefault(line[0], []).append(line[1])
    print("Read {:d} terms.".format(len(gene_dict)))
    return gene_dict


def read_results(res_file):
    results = []
    with open(res_file, 'r') as file:
        res_reader = csv.reader(file)
        for line in res_reader:
            results.append(line)
    return results


def translate_results(results, gene_dict, gene_length):
    for result in results:
        raw_pattern = result[0].strip()
        genes = [raw_pattern[i:i+gene_length] for i in range(0, len(raw_pattern), gene_length)]        
        [print(gene_dict[gene][0], end='') for gene in genes]
        print(',',','.join(result[1:]),sep='')


if len(sys.argv) is not 4:
    print("Usage: {} <results file> <dictionary file> <gene_length>".format(sys.argv[0]))
    exit(1)
else:
    results_file = sys.argv[1]
    dict_file = sys.argv[2]
    gene_length = int(sys.argv[3])

gene_dict = read_dict(dict_file)
results = read_results(results_file)
translation = translate_results (results, gene_dict, gene_length)
