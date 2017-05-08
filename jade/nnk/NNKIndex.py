import sys, re

import jade.basic.path
from collections import defaultdict, namedtuple

from pyigclassify.IgClassifyFASTA import IgClassifyFASTA
from pyigclassify.IgClassifyPDB import IgClassifyPDB
from pyigclassify.modules.chains.util import ClassifiedAb
from pyigclassify.modules.chains.AbChain import AbChain
from pyigclassify.modules.Structure import PDBResInfo
from pyigclassify.tools.path import open_file
from pyigclassify.modules.chains.IgChain import light_genes
from jade.nnk.NNKAbMaturation import GetNNKData

#Chaintypes
_heavy_ = 0
_light_ = 1

_chain_id_ = defaultdict()
_chain_id_[_heavy_] = 'heavy'
_chain_id_[_light_] = 'light'


#Skip CDRs: List of CDRs
#Skip Residues: List of Residue tuples to skip (pdb_res, chain, icode(' '))
#Skip Residues Number: List of indexes.
NumericalIndexOptions = namedtuple('NumericalIndexOptions', 'skip_cdrs skip_residues_tuple skip_residues_number include_only_regions')

class TemplateAbNNKIndex(object):
    """
    Reads the relative enrichment CSV file from our NNK data, parses out the sequence from this,
     and holds PyIgClassify identified chains types

    """
    def __init__(self, class_ab = "glCHA31", data_type = _heavy_, options = None):

        self.chainID = _chain_id_[data_type]
        self.class_ab_ = class_ab
        self.chain, self.res_info= self.__setup_from_class_directory(data_type)

        if not options:
            opt = NumericalIndexOptions(['H3'], [], [],['FRAME', 'CDR']) #Default indexing options
            self.__set_numerical_index_options(opt)
        else:
            self.__set_numerical_index_options(options)

        self.index, self.reverse_index = self.__setup_indexing()
        #print repr(self.index)

    def __set_numerical_index_options(self, opt):
        """
        Set a namedtuple with specific options for indexing.
        :type options: NumericalIndexOptions
        :return:
        """
        self.index_options = opt
        #self.index, self.reverse_index = self.__setup_indexing()


    def __setup_from_class_directory(self, chain_type=_heavy_):
        """
        Sets up the template index from a germline sequences file in the database.

        :param chain_type:
        :rtype: ClassifiedAb, PDBResInfo
        """
        germline_directory = defaultdict()
        inpath = jade.basic.path.get_nnk_database_path()
        print "Reading "+inpath+"/sort_sequences.txt"
        INFILE = open(inpath + "/sort_sequences.txt", 'r')
        for line in INFILE:
            line = line.strip()
            if not line or line[0] == '#': continue

            lineSP = line.split()
            #print line
            germline_directory[lineSP[0]] = lineSP[1]
        INFILE.close()

        fasta_path = inpath+"/"+germline_directory[self.class_ab_]
        classifier = IgClassifyFASTA(fasta_path)
        chains = classifier.run(write_data=True)
        chain = chains[0]
        if not isinstance(chain, ClassifiedAb): sys.exit()

        ab_index = chain.ab_chain.get_pdb_res_info()

        print "Index Information:"
        print str(chain.ig_chain)
        #print chain.ab_chain.get_fasta_print()

        return chain, ab_index

    def __setup_indexing(self):
        index = defaultdict()
        if not isinstance(self.res_info, PDBResInfo): sys.exit()

        print repr(self.res_info.total_residues())
        for i in range(1, self.res_info.total_residues() +1 ):
            res = self.res_info.get_residue(i)
            tup = res.tuple()

            if self.index_options.skip_cdrs and res.is_cdr() and res.get_cdr_type() in self.index_options.skip_cdrs:
                print "Skipping "+res.get_cdr_type()+" From Indexing"
                continue
            elif self.index_options.skip_residues_tuple and tup in self.index_options.skip_residues_tuple:
                print "Skipping set skip: "+repr(tup)+" From Indexing"
                continue
            elif self.index_options.skip_residues_number and i in self.index_options.skip_residues_number:
                print "Skipping set resnum: "+repr(tup)+" From Indexing"
                continue
            elif self.index_options.include_only_regions and res.get_region() not in self.index_options.include_only_regions:
                print "Skipping Non set region "+repr(tup)
                continue

            index[tup] = i

        reverse_index = defaultdict()
        for key, value in index.iteritems():
            reverse_index[value] = key

        return index, reverse_index


class TestAbNNKIndex(object):
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







