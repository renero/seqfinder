import logging
from seqlogging import setup_logging
import getopt
import sys
import os

from core import best_distinctive
from seqhandling import show_dss, save_dss
from seqread import read_seq
import settings


def checkargs(argv):
    outputfile_path = ''
    usage = '{} [-h][-o outputfile]'.format(os.path.basename(__file__))
    try:
        opts, args = getopt.getopt(argv, "ho:", ["ofile="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit(2)
        elif opt in ("-o", "--ofile"):
            outputfile_path = arg
            if os.path.isfile(outputfile_path):
                logger.error("Output file already exists. Use another filename.")
                sys.exit(3)

    return outputfile_path


settings.init()
setup_logging()
logger = logging.getLogger(__name__)
outputfile = checkargs(sys.argv[1:])

data_bajas = read_seq(settings.bajas, settings.delimiter)
data_altas = read_seq(settings.altas, settings.delimiter)

dss = []
logger.debug("Number of iterations = {:d}".format(settings.NOI))
logger.debug("Number of samples = {:d}".format(settings.NOS))
logger.debug("Confidence level = {:f}".format(settings.MCP))
logger.debug("Max nr. of wildcards accepted = {:d}".format(settings.HMW))
logger.debug("Exhaustive search = {:s}".format(str(settings.exhaustive_search)))

for los in range(settings.min_seq_len, settings.max_seq_len, settings.gene_length):
    logger.debug("Current Length-of-sequence = {:d}".format(los))
    best_distinctive(dss, data_bajas, data_altas,
                     num_samples=settings.NOS,
                     ntimes=settings.NOI,
                     seq_length=los,
                     regex_length=settings.HMW,
                     seq_min_count_percentage=settings.MCP,
                     min_g_score=settings.min_g_score,
                     filter_distinctive=settings.filter_distinctive,
                     exhaustive_search=settings.exhaustive_search,
                     gene_length=settings.gene_length)
if len(dss) is 0:
    logger.warning('List of sequences appears to be empty.')
else:
    show_dss(dss)
    if len(sys.argv) > 1:
        logger.debug("Saving output to CSV file ({})".format(outputfile))
        save_dss(dss, outputfile)
