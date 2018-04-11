import pandas as pd


def read_seq(filename=None, delimiter=','):
    """
    Read a text CSV filename and returns the contents stored in-memory.
    :param filename: the name of the file
    :return: the contents of the second column in the CSV.
    """
    if filename is None:
        filename = "/Users/jrenero/Code/distinctive_sequences/sample.txt"
    data = pd.read_csv(filename, header='infer', sep=delimiter)
    return data.iloc[:, 1].values.tolist()
