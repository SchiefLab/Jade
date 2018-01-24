import os,sys,re
from collections import defaultdict
from jade.basic.path import open_file
from jade.basic.path import get_decoy_name

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
            lineSP = line.split()
            if len(lineSP) == 1:
                data['decoy'] = get_decoy_name( os.path.basename(decoy_path) )
            else:
                data['decoy'] = get_decoy_name( lineSP[-1] )

        elif line.startswith("label"):
            labels = line.split()[1:]
            labels = ["total_score" if x=="total" else x for x in labels ]

        elif line.startswith("pose"):
            scores = [ float(x) for x in line.split()[1:] ]

        else:
            continue

    INFILE.close()

    for i in range(0, len(labels)):
        data[labels[i]] = scores[i]

    return data