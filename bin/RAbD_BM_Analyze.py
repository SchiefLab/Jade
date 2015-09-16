#!/usr/bin/env python

from Tkinter import *
#from ttk import *

import tkFileDialog
import tkMessageBox
import tkSimpleDialog

import os
import sys
import re
import copy
from collections import defaultdict
from argparse import ArgumentParser

import datetime

from RAbD_BM.BenchmarkInfo import *
from tools.general import *
from tools.path import *

import RAbD_BM.PoolData as pool
import RAbD_BM.tools as tools

from bin import create_features_json as create_json
from glob import glob

##This was definitely not the best way to do this, but its here and its too late to reimagine.


class AnalyzeBenchmarks:
    def __init__(self):
        """
        full_name: $exp.$pdbs.$type.$benchmark
        final_name: docked_relax_top
        """
        self.base_dir = os.getcwd()
        self._parse_options()

        self.native_db = self.base_dir+"/"+self.options.bm_set_db_dir+"/"+"natives."+self.options.bm_set+"."+ "all"+"."+self.options.scorefunction+".db3"
        self.lambda_db = self.base_dir+"/"+self.options.bm_set_db_dir+"/"+"natives."+self.options.bm_set+"."+ "lambda"+"."+self.options.scorefunction+".db3"
        self.kappa_db = self.base_dir+"/"+self.options.bm_set_db_dir+"/"+"natives."+self.options.bm_set+"."+ "kappa"+"."+self.options.scorefunction+".db3"

        self.recovery_base = self.base_dir+"/features/recovery"

        if not self.options.paper_ab_db:
            self.ab_db = os.getenv("ROSETTA3_DB")+"/sampling/antibodies/antibody_database_rosetta_design.db"
        else:
            self.ab_db = os.getenv("ROSETTA3_DB") +"/sampling/antibodies/antibody_database_rosetta_design_north_paper.db"

        if not os.path.exists(self.native_db):
            sys.exit("Native DB path does not exist: "+self.native_db)

        self.date = tools.get_today()
        self.extension = get_rosetta_program("", False, self.options.compiler)

        self._setup_benchmark_from_options()

    def _parse_options(self):
        ####################################################################################################################
        ##                                                  OPTIONS
        ####################################################################################################################
    
    
        parser = ArgumentParser()
    
        ############################
        ## Required Options
        ############################
        parser.add_argument("--decoy_dir",
                          help = "rel path to where decoys exist for benchmark")
    
    
        parser.add_argument("--skip_antibody_features",
                          help = "Skip output of the AntibodyFeatures db (Recoveries are not calculated for the AntibodyFeatures runs.",
                          default = False,
                          action = "store_true")
    
        parser.add_argument("--skip_cluster_features",
                          help = "Skip output of the cluster features db (which includes calculating recoveries)",
                          default = False,
                          action = "store_true")
    
        parser.add_argument("--paper_ab_db",
                          help = "Used paper ab db for analysis",
                          default = True)
    
        #######################
        ## DB input options
        #######################
    
        parser.add_argument("--bm_set",
                          default="new_20",
                          help = "The benchmark set of the input databases.  Will use this native to fine the correct input databases for recovery comparisons.")
    
        parser.add_argument("--bm_set_db_dir",
                          default="databases/natives",
                          help = "Path to the set of native databases for recovery.")
    
    
        ############################
        ## Single experiment options
        ############################
        parser.add_argument("--full_name",
                          help = "Full name of the experiment: $exp.$pdbs.$type.$benchmark")
    
        parser.add_argument("--final_name",
                          help = "Final name used to report data: docked_relax_top")
    
    
    
    
        #######################
        ## Not Required options
        #######################
        parser.add_argument("--scorefunction",
                          help = "The scorefunction used for features reporters and analysis",
                          default = "talaris2013")
    
        parser.add_argument("--exp_batch",
                          help = "Name of the batch to use for databases",
                          default = "FinalPaperBM.8.2015")
    
    
    
        parser.add_argument("--get_ensemble_data",
                          help = "Attempt to get data for an ensemble as well",
                          default = False,
                          action = "store_true")
    
        parser.add_argument("--check_nstruct",
                          help = "Make sure all decoys are present",
                          default = False,
                          action = "store_true")
    
        parser.add_argument("--delete_dbs",
                          help = "Delete output databases if already present instead of skipping them.",
                          default = False,
                          action = "store_true")
    
        parser.add_argument("--compiler",
                          default = "clang")

        self.options = parser.parse_args();

    def _setup_benchmark_from_options(self):
        BM_info = BenchmarkInfo()
        BM_info.set_data(self.options.full_name,
                         self.options.final_name,
                         self.options.scorefunction,
                         self.options.decoy_dir,
                         self.options.exp_batch)

        self.set_benchmarks([BM_info])



    def _setup_benchmarks(self):
        """
        Checks benchmarks for expected and found decoys.  Returns ones that pass cutoffs if we are checking nstruct.
        """
        final_benchmarks = []
        for benchmark in self.benchmarks:
            if isinstance(benchmark, BenchmarkInfo):pass

            if self._benchmark_ok(benchmark):
                final_benchmarks.append(benchmark)
            else:
                print "Something is wrong with benchmark.  Either not enough structures or paths do not match!  Skipping..."
                print benchmark.final_name+" "+benchmark.full_name
                continue

        return final_benchmarks

    def _benchmark_ok(self, benchmark):
        """
        Checks the total number of decoys and paths, writes to analysis log.  Returns boolean for success
        """

        if not os.path.exists(benchmark.decoy_dir):
            print "Benchmark decoy dir not found!"
            print benchmark.decoy_dir
            return False

        else:
            return True

    def add_benchmark(self, benchmark):
        if isinstance(benchmark, BenchmarkInfo): pass
        self.benchmarks.append(benchmark)

    def set_benchmarks(self, benchmarks):
        self.benchmarks = benchmarks


    def set_delete_dbs(self, delete_dbs):
        """
        Set a boolean to delete the present dbs.
        If True, we delete and rerun.
        If False, we keep and skip.
        """
        self.options.delete_dbs = delete_dbs

    def analyze_benchmarks(self):
        """
        Runs the FUll analysis
        """
        benchmarks = self._setup_benchmarks()

        for benchmark in benchmarks:
            self._run_analysis(benchmark, self.options.get_ensemble_data, self.options.skip_antibody_features, self.options.skip_cluster_features)

    def _run_analysis(self, benchmark, get_ensemble_data = False, skip_antibody_features = False, skip_cluster_features = False):
        """
        1) Makes a PDBList of the decoys.  Runs features reporter to collect recovery data
        2) Creates the JSON to run the Features R scripts
        3) Analyzes the features database for recovery by running pool data.
        """

        def run_features(features_type):
            db_path = self.base_dir+"/databases/"+benchmark.full_name+".top."+features_type+"."+benchmark.scorefunction+".db3"
            print db_path
            #os.remove(db_path)
            features_cmd = "./make_pdblists_run_features.sh "+benchmark.decoy_dir+" "+benchmark.full_name+" "+features_type+" "+benchmark.scorefunction+" "+benchmark.exp_batch + " "+repr(get_ensemble_data)+" "+self.extension

            if not os.path.exists(db_path):
                print "Making PDBLists and running Features Reporters for "+features_type
                print features_cmd
                os.system(features_cmd)

            elif self.options.delete_dbs:
                os.remove(db_path)

                print "Making PDBLists and running Features Reporters for "+features_type
                print features_cmd
                os.system(features_cmd)
            else:
                print db_path+" exists...skipping"

            if features_type == "cluster_features":
                self._calculate_recoveries(db_path, benchmark)
                self._pool_data(db_path, benchmark)


                if get_ensemble_data:
                    print "Calculating ensemble data"
                    db_path = self.base_dir+"/databases/"+benchmark.full_name+".ens.cluster_features."+benchmark.scorefunction

                    if os.path.exists(db_path):
                        os.remove(db_path)

                    new_benchmark = copy.deepcopy(benchmark)
                    new_benchmark.final_name = benchmark.final_name+".ens"

                    self._calculate_recoveries(db_path, benchmark)
                    self._pool_data(db_path, benchmark, False)




        if isinstance(benchmark, BenchmarkInfo): pass

        add_rec_data_to_current_db = True

        if not skip_cluster_features:
            run_features("cluster_features")

        if not skip_antibody_features:
            run_features("antibody_features")

    def _calculate_recoveries(self, db_path, benchmark):
        """
        Creates JSON, runs R scripts to calculate recoveries.  Moves tables to final output directory
        """

        print "Making jsons and running recovery R script..."

        if os.path.exists(self.base_dir+"/build"):
            os.system("rm -r "+self.base_dir+"/build")

        create_json.write_json_for_single_recovery_experiment(db_path, self.native_db, benchmark.final_name, self.base_dir+"/features")
        json_file = self.base_dir+"/features/jsons/cluster_features."+benchmark.final_name+".json"

        r_cmd = get_rosetta_features_root()+"/compare_sample_sources.R --config "+json_file
        print r_cmd

        new_recovery_path = self.recovery_base+"/"+benchmark.final_name
        if os.path.exists(new_recovery_path):
            os.system("rm -r "+new_recovery_path)

        os.system(r_cmd)
        os.system("mkdir "+new_recovery_path)
        #os.system(r_cmd)

        #Move cdr cluster features since I still can't seem to fix giving an output dir and output name for R features reports.
        #This makes it so running in parellel is difficult since we may be overwriting each other!!!!!

        print "Moving recovery tables to final directory"

        dirs = ["output_csv", "output_html", "output_tab_delimited_table"]
        for d in dirs:
            print("cp -r"+self.base_dir+"/build/cdr_cluster_recovery/"+d+" "+new_recovery_path)
            os.system("cp -r "+self.base_dir+"/build/cdr_cluster_recovery/"+d+" "+new_recovery_path)

        os.system("rm -r "+self.base_dir+"/build")

    def _pool_data(self, db_path, benchmark, top = True):
        """
        Runs Pool data
        """
        print "Pooling Data"

        if top:
            recovery_db_name = benchmark.full_name+".top.recoveries."+benchmark.scorefunction+".db3"

        pool_data = pool.PoolData(self.lambda_db, self.kappa_db, self.ab_db, self.base_dir+"/features", [benchmark])
        pool_data.apply(self.base_dir+"/databases", recovery_db_name,  True)

        #Add to current db:
        pool_data.output_data(True, os.path.dirname(db_path), os.path.basename(db_path), True) #Drop tables if exist.

    ####################################################################################################################

    def _write_exp_log(self, benchmark, log_file = "benchmark_list.txt"):
        """
        Writes a simple log file that will be used later to control plotting and data analysis.
        """
        if not os.path.exists(log_file):
            FILE = open(log_file, 'w')
            FILE.write("#date full_name final_name scorefunction expected_nstruct found_nstruct exp_batch rel_decoy_path rel_log_path\n")
        else:
            FILE = open(log_file, 'a')

        line = (benchmark.date+"\t"+benchmark.full_name+"\t"+benchmark.final_name+"\t"+benchmark.scorefunction+"\t"+repr(benchmark.expected_nstruct)+"\t"+repr(benchmark.found_nstruct)+"\t"+
                   benchmark.exp_batch+"\t"+benchmark.decoy_dir+"\t"+benchmark.log_dir+"\n")

        #print(line)
        FILE.write(line)

        FILE.close()


#############################################
### Controls ALL analysis of decoys and data.
#############################################
if __name__=="__main__":

    analyzer = AnalyzeBenchmarks()
    analyzer.analyze_benchmarks()

    print "Complete"