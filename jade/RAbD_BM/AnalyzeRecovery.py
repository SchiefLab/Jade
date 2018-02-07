#!/usr/bin/python
#Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#This Script parses the log of MPI-run antibody_design benchmarks to get the frequency of clusters and lengths against a native database.
#It also gets DB frequences from our main databases and will soon parse the resultant length and cluster recoveries to determine the 
# Risk Ratio of each experiment for conclusive statistics.


import gzip,os,sys,re,glob,pandas


try:
    import sqlite3
    import numpy
except ImportError:
    print "Cannot import SQLITE3 or NUMPY.  Pleas install to use pool data"

from collections import defaultdict
from optparse import OptionParser

import jade.RAbD_BM.tools_ab_db as pyig_tools
import jade.RAbD_BM.tools_features_db as feat_tools
from jade.rosetta_jade.BenchmarkInfo import BenchmarkInfo
from jade.RAbD_BM.AnalysisInfo import AnalysisInfo
from jade.RAbD_BM.AnalysisInfo import NativeInfo
from jade.basic.structure.Structure import AntibodyStructure
from jade.basic.structure.Structure import cdr_names
from jade.basic.structure.Structure import heavy_cdr_names

from jade.basic.numeric import *

class AnalyzeRecovery:
    """
    Pools Recovery and RR data, outputs to DB
    """
    def __init__(self, pyig_design_db_path, analysis_info, native_info, cdrs = None):
        """
        Class for anlyzing length and cluster recovery and risk ratios.

        Optionally pass in CDRs to control which CDRs are analyzed.
          Otherwise, we read the info from RUN_SETTINGS.txt (which is used for benchmarking), which either has ALL or a specific CDR.

        
        :param pyig_design_dbpath: str
        :param analysis_info: AnalysisInfo
        :param native_info: NativeInfo
        :param cdrs: [str]
        """

        cdrs = ["L1", "L2", "L3", "H1", "H2"]


        if not isinstance(analysis_info, AnalysisInfo): sys.exit()
        if not isinstance(native_info, NativeInfo): sys.exit()

        self.pyig_db_path = pyig_design_db_path
        self.analysis_info = analysis_info
        self.native_info = native_info
        self.lambda_kappa_dict = defaultdict()


        self.lambda_kappa_dict["lambda"] = native_info.lambda_pdbids
        self.lambda_kappa_dict["kappa"] = native_info.kappa_pdbids

        ## Set the CDRs we will need. (WE NEED TO MAKE SURE WE HAVE H# OR NOT!!!)
        if cdrs:
            self.cdrs = cdrs
        elif self.analysis_info.bm_info.settings["CDR"] == "ALL":
            self.cdrs = cdr_names
        else:
            self.cdrs = [self.analysis_info.bm_info.settings["CDR"]]


        self.recovery_types = ["length", "cluster"]

        ##Create empty result dataframes.
        #Final Data
        self.recovery_df = pandas.DataFrame()
        self.observed_df = pandas.DataFrame()
        self.represented_df = pandas.DataFrame()

        ##Initialize the calculators
        self.recovery_calc = TopRecoveryCalculator(self.native_info.get_features_db())
        self.observed_calc = ObservedRecoveryCalculator(self.native_info.get_features_db())

        self.initialize()

    def initialize(self):
        """
        Initialize ALL input data before calculating and outputing everything.
        """
        if not os.path.exists(self.native_info.get_features_db()):
            sys.exit("native_db path not good: "+self.native_info.get_features_db())


        self.recovery_df = self.recovery_calc.apply(
            self.analysis_info.get_exp(),
            self.native_info.pdbids,
            self.cdrs,
            self.analysis_info.get_features_db())


        self.observed_df = self.observed_calc.apply(
            self.analysis_info.get_exp(),
            self.native_info.pdbids,
            self.cdrs,
            self.analysis_info.get_decoy_dir())

    def apply(self, db_path, drop_tables = False):
        """
        Calculate and Output all the data to the given database.

        :param db_path: str
        """
        print "opening "+db_path
        db_con = sqlite3.connect(db_path)


        table_action = "replace" if drop_tables else "append"

        if drop_tables:
            for table in ["full_data", "by_cdr_data", "by_exp_data"]:
                db_con.execute("DROP TABLE IF EXISTS "+table)


        full_df = calculate_recovery_and_risk_ratios(self.recovery_df,self.observed_df)
        cdr_df = calculate_per_cdr_rr_and_recovery(self.analysis_info.get_exp(), self.cdrs, full_df)
        exp_df = calculate_exp_rr_and_recovery(self.analysis_info.get_exp(), full_df)

        #Write them all out to sql
        full_df.to_sql("full_data", db_con, if_exists=table_action)
        cdr_df.to_sql("by_cdr_data", db_con, if_exists=table_action)
        exp_df.to_sql("by_exp_data", db_con, if_exists=table_action)

        db_con.close()

######### Recovery Calculators

class RecoveryCalculator(object):
    def __init__(self, native_db_path):
        """
        Super simple base class for recovery calculators.

        """
        self.all_cdrs = cdr_names
        self.heavy_cdrs = heavy_cdr_names
        self.Antibody = AntibodyStructure()
        print "Opening "+native_db_path
        self.native_df = feat_tools.get_cdr_cluster_df(native_db_path)
        #print "Native Dataframe: "
        #print self.native_df

    #IMPLEMENT THIS:
    #def apply(self):
    #    """
    #    Each calculator implements its own version of apply (which may take different arguments.
    #      Pythons argument dictionary is horrible, and I hate it, so instead, this is a comment
    #
    #    Each calculator also returns the resultant dataframe.
    #
    #    :return: pandas.DataFrame
    #    """
    #    pass

class TopRecoveryCalculator(RecoveryCalculator):
    def __init__(self, native_db_path):
        """
        Calculate length and cluster recovery from the features database and native database,
         which has the cluster and length info.
        """
        RecoveryCalculator.__init__(self, native_db_path)

        self.bm_df = pandas.DataFrame

    def apply(self, exp_name, pdbids, cdrs, bm_db_path, output_dir = "data" ):
        """
        Calculate length and cluster recoveries.  Store them the same way we used to for the recovery parser.
        Returns the resulting dataframe of recoveries.
        :rtype: pandas.DataFrame

        """

        #Reset Data.
        self.bm_df = pandas.DataFrame

        #Calculate
        exp_name = exp_name

        print "Getting BM data: "+bm_db_path

        if not os.path.exists(bm_db_path):
            sys.exit("BM DB does not exist!")

        bm_df = feat_tools.get_cdr_cluster_df(bm_db_path)

        #print self.native_df
        flat_dict = defaultdict(list)
        for pdbid in pdbids:
            for cdr in cdrs:
                print pdbid+" "+cdr

                native_length = feat_tools.get_length(self.native_df, pdbid, cdr)
                native_cluster = feat_tools.get_cluster(self.native_df, pdbid, cdr)

                #print "Native: "+repr(native_length)
                #print "Cluster: "+repr(native_cluster)

                total_entries = feat_tools.get_total_entries(bm_df, pdbid, cdr)

                length_recovery = feat_tools.get_length_recovery(bm_df, pdbid, cdr, native_length)
                cluster_recovery = feat_tools.get_cluster_recovery(bm_df, pdbid, cdr, native_cluster)

                #print "length_recovery: "+str(length_recovery)
                #print "cluster_recovery: "+str(cluster_recovery)

                flat_dict['length_recovery_freq'].append(length_recovery)
                flat_dict['cluster_recovery_freq'].append(cluster_recovery)
                flat_dict['total_entries'].append(total_entries)
                flat_dict['pdbid'].append(pdbid)
                flat_dict['cdr'].append(cdr)
                flat_dict['exp'].append(exp_name)


        self.bm_df = bm_df
        result_df = pandas.DataFrame.from_dict(flat_dict)

        ##Convert into something easily converted to pandas and output CSV.

        print "Recoveries calculated."

        self.native_df.to_csv(output_dir+"/native_clusters.csv")

        columns = ['exp', 'pdbid', 'cdr', 'length_recovery_freq', 'cluster_recovery_freq', 'total_entries']

        out = output_dir+"/all_recoveries.csv"
        if os.path.exists(out):
            ftype = 'a'
            header = False
        else:
            ftype = 'w'
            header = True

        OUTFILE = open(out, ftype)
        result_df.to_csv(OUTFILE, columns=columns, header=header)
        OUTFILE.close()

        return result_df

class ObservedRecoveryCalculator(RecoveryCalculator):
    def __init__(self, native_db_path):
        """
        Calculates the number of times the native clusters and lengths were observed during the experiment, for each PDB.
        """
        RecoveryCalculator.__init__(self, native_db_path)

    def apply(self, exp_name, pdbids, cdrs, bm_decoy_path, output_dir = "data"):
        """
        Calculates the number of times the native clusters and lengths were observed during the experiment, for each PDB.
        Returns the resulting dataframe.

        :rtype: pandas.DataFrame
        """

        flat_dict = defaultdict(list)
        data_types = ["total_grafts", "native_lengths_observed", "native_clusters_observed"]
        for pdbid in pdbids:
            print pdbid
            totals = defaultdict()
            clusters = defaultdict()
            lengths = defaultdict()

            ##Initialize native lengths and clusters to do the check.
            for cdr in cdrs:

                clusters[cdr] = feat_tools.get_cluster(self.native_df, pdbid, cdr)
                lengths[cdr] = feat_tools.get_length(self.native_df, pdbid, cdr)

                ##Initialize the dictionary that we will use to do the counts.  This dictionary is nessessary to do the counts!
                totals[cdr] = defaultdict()
                totals[cdr]["total_grafts"] = 0
                totals[cdr]["native_lengths_observed"] = 0
                totals[cdr]["native_clusters_observed"] = 0

            filenames = get_decoys(bm_decoy_path, pdbid)
            for filename in filenames:
                total_grafts = 0
                if os.path.basename(filename).split('.')[-1] == "gz":
                    INFILE = gzip.open(filename, 'rb')
                else:
                    INFILE = open(filename, 'r')

                for line in INFILE:
                    line = line.strip()
                    if not line or line.startswith("ATOM"): continue

                    if re.search(" DATA GRAFT_CLOSURE ", line):
                        total_grafts+=1
                        lineSP = line.split()
                        data_index = 0
                        for i in range(0, len(lineSP)):
                            if lineSP[i]=="DATA":
                                data_index = i
                                break

                        input_tag = lineSP[data_index+3]
                        cluster = lineSP[data_index+5]
                        cdr = cluster.split("-")[0]
                        length = int(cluster.split("-")[1])


                        # Now, Add the data to self.exp_data
                        #print "Adding data for input_tag "+input_tag
                        #print "Adding data for cdr " +cdr

                        totals[cdr]["total_grafts"] += 1

                        if length == lengths[cdr]:
                            totals[cdr]["native_lengths_observed"] += 1

                        if cluster == clusters[cdr]:
                            totals[cdr]["native_clusters_observed"] += 1

                #print "grafts: "+repr(total_grafts)
                INFILE.close()

            ###Add to flattened dic.
            for cdr in cdrs:

                flat_dict["pdbid"].append(pdbid)
                flat_dict["cdr"].append(cdr)
                flat_dict["exp"].append(exp_name)

                for data_type in data_types :
                    #print pdbid
                    #print cdr
                    #print exp_name

                    flat_dict[data_type].append(totals[cdr][data_type])


        ##Now we unflatten dictionary and turn it into a nice dataframe.
        columns = ["exp", "pdbid", "cdr"]
        columns.extend(data_types)

        result_df = pandas.DataFrame.from_dict(flat_dict)
        out = output_dir+"/all_observed.csv"
        if os.path.exists(out):
            ftype = 'a'
            header = False
        else:
            ftype = 'w'
            header = True

        OUTFILE = open(out, ftype)
        result_df.to_csv(OUTFILE, columns=columns, header=header)
        return result_df

class PyIgClassifyDBRepresentationCalculator(RecoveryCalculator):
    def __init__(self, native_db_path):
        """
        Calculates the number of times lengths and clusters are present in the PyIgClassify database.
        :return:
        """
        RecoveryCalculator.__init__(self, native_db_path)

    def apply(self, exp_name, cdrs, pyig_db_path, lambda_kappa_dict, output_dir = "data"):
        """
        Calculates the number of times lengths and clusters are present in the PyIgClassify database.
        :param lambda_kappa_dict : dict-like ['lambda'] = [pdbid,]

        :rtype: pandas.DataFrame
        """
        flat_dict = defaultdict(list)
        cdr_data_df = pyig_tools.get_cdr_data_table_df(pyig_db_path)
        for light_gene in lambda_kappa_dict:
            for pdbid in lambda_kappa_dict[light_gene]:
                for cdr in cdrs:

                    gene = "heavy" if cdr in self.heavy_cdrs else light_gene
                    native_length = feat_tools.get_length(self.native_df, pdbid, cdr)
                    native_cluster = feat_tools.get_cluster(self.native_df, pdbid, cdr)

                    length_enrichment = pyig_tools.get_length_enrichment(cdr_data_df, gene, cdr, native_length)
                    cluster_enrichment = pyig_tools.get_cluster_enrichment(cdr_data_df, gene, cdr, native_cluster)

                    total_entries = pyig_tools.get_total_entries(cdr_data_df, gene, cdr)

                    flat_dict["exp"].append(exp_name)
                    flat_dict["pdbid"].append(pdbid)
                    flat_dict["cdr"].append(cdr)
                    flat_dict["gene"].append(gene)
                    flat_dict["length_representation"].append(length_enrichment)
                    flat_dict["cluster_representation"].append(cluster_enrichment)
                    flat_dict["total_entries_of_gene_and_cdr"].append(total_entries)

        ##Now we unflatten dictionary and turn it into a nice dataframe.
        columns = ["exp", "pdbid", "gene", "cdr", "length_representation", "cluster_representation", "total_entries_of_gene_and_cdr"]
        result_df = pandas.DataFrame.from_dict(flat_dict)
        out = output_dir+"/pyigclassify_db_representation.txt"
        if os.path.exists(out):
            ftype = 'a'
            header = False
        else:
            ftype = 'w'
            header = True

        OUTFILE = open(out, ftype)
        result_df.to_csv(OUTFILE, columns=columns, header=header)
        return result_df


####### Combine data and get risk ratios, recovery percentage, etc.

def calculate_recovery_and_risk_ratios(top_recovery_df, observed_df):
    """
    Calculate the Risk Ratio and Recovery Percent for each pdb/cdr given dataframes output by the calculators below.

    Return a merged dataframe of the top recovery and observed, with the resulting risk ratio data.

    :param top_recovery_df: pandas.DataFrame
    :param observed_df: pandas.DataFrame
    :rtype: pandas.DataFrame
    """
    df = pandas.merge(top_recovery_df, observed_df, how="outer")

    for rtype in ["length", "cluster"]:

        df[rtype+'_recovery'] = (df[rtype+'_recovery_freq']/df['total_entries'].astype('float'))
        df[rtype+'_rr'] = df[rtype+'_recovery']/(df['native_'+rtype+'s_observed']/df['total_grafts'].astype('float'))

    return df

def calculate_per_cdr_rr_and_recovery(exp, cdrs, result_df):
    """
    Calculate the recovery and risk-ratios PER CDR.
    :rtype: pandas.DataFrame
    """
    flat_dict = defaultdict(list)
    for cdr in cdrs:
        flat_dict['exp'].append(exp)
        flat_dict['cdr'].append(cdr)
        for rtype in ['length', 'cluster']:

            flat_dict[rtype+'_rr'].append(numpy.mean(result_df[(result_df['cdr'] == cdr)][rtype+'_rr']))
            flat_dict[rtype+'_recovery'].append(numpy.mean(result_df[(result_df['cdr'] == cdr)][rtype+'_recovery']))

    return pandas.DataFrame.from_dict(flat_dict)

def calculate_exp_rr_and_recovery(exp, result_df):
    """
    Calculate the overall recovery and risk ratio.
    :param exp:
    :param result_df:
    :rtype: pandas.DataFrame
    """
    flat_dict = defaultdict(list)
    flat_dict["exp"].append(exp)
    for rtype in ['length', 'cluster']:
        flat_dict[rtype+'_rr'].append(numpy.mean(result_df[result_df['exp'] == exp][rtype+'_rr']))
        flat_dict[rtype+'_recovery'].append(numpy.mean(result_df[result_df['exp'] == exp][rtype+'_recovery']))

    return pandas.DataFrame.from_dict(flat_dict)

#########

def get_decoys(input_dir, pdbid):
    """
    Use GLOB to Match on pdbid for file names in the input dir.
    This should skip all the extra PDBs like excn, initial, relax, etc.
    :param input_dir: str
    :param tag: str
    """
    if not os.path.exists(input_dir):
        sys.exit("MPI Log Dir does not exist!")
    search_name = input_dir+"/*"+pdbid+"*.pdb*"
    #print search_name

    final_filenames = []
    names = glob.glob(search_name)
    for pdb in names:
        if re.search("initial_benchmark_perturbation", pdb) or re.search("excn", pdb):
            continue
        else:
            final_filenames.append(pdb)

    return final_filenames







