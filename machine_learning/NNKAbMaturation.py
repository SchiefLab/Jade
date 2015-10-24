#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas
import glob
from collections import defaultdict

import machine_learning.util as ml_util



sorts=['S1','S2','S3']

class GetNNKData(object):
    """
    Get NNK Data as a formatted tupple of 1d data (Or raw Pandas DF)
    """
    def __init__(self, data_dir, ab_class="glCHA31" ):
        self.data_dir = data_dir
        self.class_dir = data_dir+"/"+ab_class
        self.ab_class = ab_class
        self.data_types = self.__get_all_data_types()
        self.antigens = self.__get_all_antigens()

    def reinit(self, ab_class):
        self.__init__(self.data_dir, ab_class)

    def __get_unique_data_groups(self, split_index):
        types = defaultdict()
        files = glob.glob(self.class_dir+"/*.csv")

        for f in files:
            fSP = os.path.basename(f).split("_")
            types[fSP[split_index]] = 0
        return [x for x in types]

    def __get_all_antigens(self):
        return self.__get_unique_data_groups(data_dir, 2)

    def __get_all_data_types(self):
        return self.__get_unique_data_groups(data_dir, 0)

    def get_nnk_data(self, dt, antigen="C5-SOSIP", sort="S1"):
        """
        Get pandas dataframe of NNK data.
        """
        if not dt in self.data_types:
            print "datatype not understood"
            print "Available datatypes: "+repr(self.data_types)

        if not sort in sorts:
            print "sort not understood"

        if not antigen in self.antigens:
            print "antigens not understood"

        filename=self.class_dir+"/"+"_".join([dt, antigen, sort])+".csv"
        print filename
        df = pandas.read_csv(filename)
        df.columns = df.columns.str.replace('Unnamed: 0','ResType')
        df = df.set_index('ResType')
        return df

    def get_1d_data_tuple_freq_nnk_data(self, antigen = "C5-SOSIP", sort = "S1"):
        top_freq = self.get_nnk_data('freqTopPerPosition', antigen=antigen, sort=sort)
        bot_freq = self.get_nnk_data('freqBotPerPosition', antigen=antigen, sort=sort)

        return (top_freq/bot_freq).as_matrix().flatten()

def main_test_1(data_dir, mat_sort="S1"):
    """
    Here, this is a test bed for SVM and simple neural networks
    No recurrent Neural nets or anything fancy.  Will have to try that next.
    :param data_dir:
    :return:
    """

    #Here, we need to organized the data

    data_types = ["glCHA31", "matCHA31", "matVRC01", "minVRC01"]

    # First, we get the groups.  Then, we can figure out what to do with them.  In the beginning, we just have each as separate groups.
    data_loader = GetNNKData(data_types[0])

    germline_data_sort1 = data_loader.get_1d_data_tuple_freq_nnk_data(antigen="GT8", sort="S1")
    germline_data_sort2 = data_loader.get_1d_data_tuple_freq_nnk_data(antigen="GT8", sort="S2")
    germline_data_sort3 = data_loader.get_1d_data_tuple_freq_nnk_data(antigen="GT8", sort="S3")

    ##Add sort 2?

    data_loader.reinit(data_types[1])
    mature_data = defaultdict()


    for antigen in data_loader.antigens:
        for sort in sorts:
            mature_data[sort][antigen] = data_loader.get_1d_data_tuple_freq_nnk_data(antigen, sort)

    #Use sort 3, core as a set, GT8 as a set, and the SOSIPs as a set.

if __name__ == "__main__":


    data_dir = "/Users/jadolfbr/Documents/projects/nnk_data"
    main_test_1(data_dir)

