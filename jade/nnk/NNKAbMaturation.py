#!/usr/bin/env python

import os, re, glob
from collections import defaultdict

import pandas

_sorts_=['S1','S2','S3']
_ab_types_ = ["glCHA31", "matCHA31", "matVRC01", "minVRC01"]

class GetNNKData(object):
    """
    Get NNK Data as a formatted tupple of 1d data (Or raw Pandas DF)
    """
    def __init__(self, data_dir, ab_group="glCHA31" ):
        self.data_dir = data_dir
        self.class_dir = data_dir+"/"+ab_group
        self.ab_group = ab_group
        self.data_types = self.__get_all_data_types()
        self.antigens = self.__get_all_antigens()
        print "Antigens: "+repr(self.antigens)

        self.df = None

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
        return self.__get_unique_data_groups(2)

    def __get_all_data_types(self):
        """
        Get all the data types within an antibody group directory, such as freqTopPerPosition.
        :rtype: []
        """
        return self.__get_unique_data_groups(0)

    def get_nnk_data(self, dt="freqTopPerPosition", antigen="C5-SOSIP", sort="S1"):
        """
        Get pandas dataframe of NNK data.
        """
        if not dt in self.data_types:
            print "dt not understood"
            print "Available datatypes: "+repr(self.data_types)

        if not sort in _sorts_:
            print "sort not understood"

        if not antigen in self.antigens:
            print "antigens not understood"

        filename=self.class_dir+"/"+"_".join([dt, self.ab_group, antigen, sort])
        filename = filename+".*.csv"
        #if not os.path.exists(filename):
        #    filename = self.class_dir+"/"+"_".join([dt, self.ab_group, antigen, sort])+".csv"

        filenames = glob.glob(self.class_dir+"/*.csv")
        for f in filenames:
            if re.search(filename, f):
                df = pandas.read_csv(f)

                print "Loading "+f
                df.columns = df.columns.str.replace('Unnamed: 0','ResType')
                df = df.set_index('ResType')
                return df

        raise IOError("Could not find filename "+filename)


    def get_1d_data_tuple_freq_nnk_data(self, antigen = "C5-SOSIP", sort = "S1"):

        if not self.df.empty:
            return self.df.as_matrix().flatten()
        else:

            top_freq = self.get_nnk_data('freqTopPerPosition', antigen=antigen, sort=sort)
            bot_freq = self.get_nnk_data('freqBotPerPosition', antigen=antigen, sort=sort)
            return (top_freq/bot_freq).as_matrix().flatten()

    def get_2D_data_freq_nnk_data(self, antigen = "C5-SOSIP", sort = "S1"):
        """
        Return a dataframe with ResType as index and resnum as columns.
        :param antigen:
        :param sort:
        :rtype: pandas.DataFrame
        """
        top_freq = self.get_nnk_data('freqTopPerPosition', antigen=antigen, sort=sort)
        bot_freq = self.get_nnk_data('freqBotPerPosition', antigen=antigen, sort=sort)

        self.df = top_freq/bot_freq
        return self.df

def load_1d_data(data_dir, data_type):
    """
    Here, this is a test bed for SVM and simple neural networks
    No recurrent Neural nets or anything fancy.  Will have to try that next.

    The 1D data is so that the residuetypes all line up in the SVM.
    :param data_dir:
    :return:
    """
    loaded_data = defaultdict(dict)

    #Germline
    data_loader = GetNNKData(data_dir, data_type)

    for antigen in data_loader.antigens:
        for sort in _sorts_:
            print "Loading mature data for: "+antigen+" : "+sort
            try:
                loaded_data[sort][antigen] = data_loader.get_1d_data_tuple_freq_nnk_data(antigen, sort)
            except IOError:
                print "Data does not exist for. "+antigen+" : "+sort
                print "Continueing"
                continue

    return loaded_data

def load_2d_data(data_dir, data_type):
    """
    Load a representation of the 2D data with res and position.

    :param data_dir:
    :return:
    """
    loaded_data = defaultdict(dict)

    #Germline
    data_loader = GetNNKData(data_dir, data_type)

    for antigen in data_loader.antigens:
        for sort in _sorts_:
            print "Loading mature data for: "+antigen+" : "+sort
            try:
                loaded_data[sort][antigen] = data_loader.get_2D_data_freq_nnk_data(antigen, sort)
            except IOError:
                print "Data does not exist for. "+antigen+" : "+sort
                print "Continueing"
                continue

    return loaded_data

def write_raw_sorts(data, outdir, outname = "raw_enrichments.csv", antigen = 'GT81'):
    """
    Write the data as columns, S1, S2, S3 for easy import into matlab.

    :param data:
    :param outdir:
    :param outname:
    :param antigen:
    :return:
    """
    out = outdir + "/"+outname
    OUTFILE = open(out, 'w')
    OUTFILE.write("S1\tS2\tS3\n")

    for i in range(0, len(data['S1'][antigen])):
        line = ""
        for sort in _sorts_:
            line = line + repr(data[sort][antigen][i]) + "\t"
        line = line.strip()
        OUTFILE.write(line + "\n")
    OUTFILE.close()
    print "Done"


if __name__ == "__main__":
    data_dir = "/Users/jadolfbr/Documents/projects/nnk_data/sort_data"

    data_loader = GetNNKData(data_dir, 'glCHA31')
    d = data_loader.get_nnk_data('relativeGateEnrichment', '.*', 'S1')
    print d.tail()
