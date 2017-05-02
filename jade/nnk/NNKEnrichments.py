import sys, os, numpy
from collections import defaultdict
from copy import deepcopy
import pandas

from jade.nnk import NNKAbMaturation


class NNKEnrichments(object):
    """
    Simple class that holds all the enrichment data for a particular class, antibody, and antigen.
    """
    def __init__(self, data_dir, zeros = -2.0, class_type = 'VRC01', antibody = 'glCHA31', antigen = 'GT81', sort = 'S1'):
        """

        :param data_dir: the directory with sort data. Each set of data + antibody should be in a separate directory (Ex: glCHA31, et.c)
        :param zeros: The number we use when enrichment of top/bottom gate is zero.  This is the log(enrichment). -2.0 corresponds to an enrichment of about .08.
        :param class_type:
        :param antibody:
        :param antigen:
        :param sort:
        """

        data_loader = NNKAbMaturation.GetNNKData(data_dir, antibody)

        self.df = data_loader.get_2D_data_freq_nnk_data(antigen=antigen, sort=sort)
        self.df = self.df.applymap(numpy.log)
        self.df = self.df.replace(numpy.NINF, float(zeros))

        if not isinstance(self.df, pandas.DataFrame): sys.exit()

        self.data_1D = data_loader.get_1d_data_tuple_freq_nnk_data(antigen=antigen, sort=sort)

    def max(self, position):
        """

        Get the maximum enrichment at a particular position, and the amino acid

        :param position:
        :return:
        """

        max_index = self.df[str(position)].idxmax()
        #print repr(max_index)
        return self.value(position, max_index), max_index

    def min(self, position):
        """

        Get the minimum enrichment at a particular position, and the amino acid

        :param position:
        :return:
        """

        min_index = self.df[str(position)].idxmin()

        return self.value(position, min_index), min_index

    def value(self, position, three_letter_code):
        """

        Get the enrichment value of a particular position and code.

        :param position:
        :param three_letter_code:
        :return:
        """
        return self.df.get_value(three_letter_code, str(position))


def combine_enrichments( list_of_nnk_enrichments):
    """
    Combine a list of nnk_enrichments to populate this one.
    @type list_of_nnk_enrichments: [NNKEnrichments]
    :rtype: NNKEnrichments
    """
    print "Len: "+repr(len(list_of_nnk_enrichments))

    if len(list_of_nnk_enrichments) == 1:
        return list_of_nnk_enrichments[0]

    new_enrich = deepcopy(list_of_nnk_enrichments[0])
    dfs_2D = []
    dfs_1D = []
    for enrich in list_of_nnk_enrichments:
        #print len(enrich.df.columns)
        dfs_2D.append(enrich.df)
        dfs_1D.append(enrich.data_1D)

    df_2D = pandas.concat(dfs_2D)
    new_enrich.df = df_2D.groupby(level=0).mean()

    new_enrich.data_1D = numpy.mean(numpy.array(dfs_1D), axis = 0)
    return new_enrich