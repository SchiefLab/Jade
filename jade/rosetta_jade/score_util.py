import os,sys,re
from collections import defaultdict
from jade.basic.path import open_file

def parse_decoy_scores(decoy_path):
    """
    Parse a score from a decoy and return a dictionary. 
    :param decoy_path: 
    :return: defaultdict
    """

    data = defaultdict()
    labels = []
    scores=[]
    INFILE = open_file(decoy_path)
    for line in INFILE:
        if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
            model = line.split()[-1]
            data['decoy'] = model

        elif line.startswith("label"):
            labels = line.split()[1:]

        elif line.startswith("pose"):
            scores = [ float(x) for x in line.split()[1:] ]

        else:
            continue

    INFILE.close()

    for i in range(0, len(labels)):
        data[labels[i]] = scores[i]

    return data