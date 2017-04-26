import sys, re

import basic.path
from collections import defaultdict

from pyigclassify.IgClassifyFASTA import IgClassifyFASTA
from pyigclassify.IgClassifyPDB import IgClassifyPDB
from pyigclassify.modules.chains.util import ClassifiedAb
from pyigclassify.modules.chains.AbChain import AbChain
from pyigclassify.modules.Structure import PDBResInfo
from pyigclassify.tools.path import open_file
from pyigclassify.modules.chains.IgChain import light_genes

#Chaintypes
_heavy_ = 0
_light_ = 1

_chain_id_ = defaultdict()
_chain_id_[_heavy_] = 'heavy'
_chain_id_[_light_] = 'light'


class GermlineNNKIndex(object):
    def __init__(self, class_ab = "VRC01", chain_type = _heavy_):

        def read_germline_directory():
            germline_directory = defaultdict()
            inpath = basic.path.get_nnk_database_path()
            INFILE = open(inpath + "/germline_sequences.txt", 'r')
            for line in INFILE:
                line = line.strip()
                if not line or line[0] == '#': continue

                lineSP = line.split()

                germline_directory[lineSP[0]] = lineSP[1]
            INFILE.close()
            return germline_directory

        def load(germlines, chain_type):
            inpath = germlines[class_ab][chain_type]
            classifier = IgClassifyFASTA(inpath)
            chains = classifier.run()
            chain = chains[0]
            if not isinstance(chain, ClassifiedAb): sys.exit()

            return chain

        self.germline_directory = read_germline_directory()

        self.chains = load(self.germline_directory, chain_type)
        if not isinstance(self.chains, ClassifiedAb): sys.exit()

        self.ab_index = self.chains.get_pdb_res_info()
        if not isinstance(self.ab_index, PDBResInfo): sys.exit()



    def get_seq_up_to_h3(self):
        """
        Get an array of sequence up to H3
        :return:
        """

        if not isinstance(self.ab_index, PDBResInfo): sys.exit()

        s = []
        for i in range(1, self.ab_index.total_residues() + 1 ):
            if self.ab_index.get_residue(i).get_cdr_type() != 'H3':
                s.append(self.ab_index.get_residue(i).get_aa())
            else:
                return s

class AntibodyNNKIndex(object):
    """
    Main class to open a test sequence or structure against a score.
    """
    def __init__(self, pdb_or_fasta_path, chain_type = _heavy_):
        self.chain_type = chain_type
        self.pdb_or_fasta_path = pdb_or_fasta_path

        if re.search('pdb', pdb_or_fasta_path):
            igclassify = IgClassifyPDB(self.pdb_or_fasta_path)
            self.classified_list = igclassify.run()

        else:
            igclassify = IgClassifyFASTA(self.pdb_or_fasta_path)
            self.classified_list = igclassify.run()

        #Prune Chains that are not the chain type!
        new_list = []
        for classified in self.classified_list:
            if classified.ig_chain:

                if classified.ig_chain.is_ScFv():
                    print classified.ig_chain.id +" Is an ScFv.  Currently, we are not working with ScFvs. Continueing"
                    continue

                if classified.ab_chain:
                    if len(classified.ig_chain.ig_domains) > 1:
                        print "Unknown use case for more than one light or heavy domain.  Contueing."
                        continue

                    gene = classified.ig_chain.ig_domains[0].get_gene()
                    if self.chain_type == _heavy_ and gene == 'heavy':
                        new_list.append(classified)
                    elif self.chain_type == _light_ and gene in light_genes:
                        new_list.append(classified)
                    else:
                        print "Chain Gene does not match set chain_type. "+gene
                        continue

        self.classified_list = new_list

    def classified_list(self):
        """
        Get a list of identified chains of the desired chaintype.

        Filters out ScFvs and chains that have multiple heavy and light.
        :rtype: [ClassifiedAb]
        """

        return self.classified_list







