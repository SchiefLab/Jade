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
from optparse import OptionParser
import datetime

from python_modules.BenchmarkInfo import *
import python_modules.PoolData as pool
import python_modules.tools as tools

from features import create_features_json as create_json
from glob import glob

p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath



def main():
    ####################################################################################################################
    ##                                                  OPTIONS
    ####################################################################################################################


    parser = OptionParser()

    ############################
    ## Required Options
    ############################
    parser.add_option("--decoy_dir",
                      help = "rel path to where decoys exist for benchmark")

    parser.add_option("--expected_nstruct",
                      help = "The expected number of final decoys")

    parser.add_option("--skip_antibody_features",
                      help = "Skip output of the AntibodyFeatures db (Recoveries are not calculated for the AntibodyFeatures runs.",
                      default = False,
                      action = "store_true")

    parser.add_option("--skip_cluster_features",
                      help = "Skip output of the cluster features db (which includes calculating recoveries)",
                      default = False,
                      action = "store_true")


    #######################
    ## DB input options
    #######################

    parser.add_option("--bm_set",
                      default="new_20",
                      help = "The benchmark set of the input databases.  Will use this native to fine the correct input databases for recovery comparisons.")

    parser.add_option("--bm_set_db_dir",
                      default="databases/natives",
                      help = "Path to the set of native databases for recovery.")


    ############################
    ## Single experiment options
    ############################
    parser.add_option("--full_name",
                      help = "Full name of the experiment: $exp.$pdbs.$type.$benchmark")

    parser.add_option("--final_name",
                      help = "Final name used to report data: docked_relax_top")




    #######################
    ## Not Required options
    #######################
    parser.add_option("--scorefunction",
                      help = "The scorefunction used for features reporters and analysis",
                      default = "talaris2013")

    parser.add_option("--exp_batch",
                      help = "Name of the batch to use for databases",
                      default = "FinalPaperBM.8.2015")



    parser.add_option("--get_ensemble_data",
                      help = "Attempt to get data for an ensemble as well",
                      default = False,
                      action = "store_true")

    parser.add_option("--check_nstruct",
                      help = "Make sure all decoys are present",
                      default = False,
                      action = "store_true")

    parser.add_option("--log_dir",
                      help = "rel path to log dir for benchmark if using old-style")

    parser.add_option("--delete_dbs",
                      help = "Delete output databases if already present instead of skipping them.",
                      default = False,
                      action = "store_true")

    parser.add_option("--rosetta_extension",
                      default = "macosclangrelease")

    parser.add_option("--nstruct_tolerance",
                      default = 5)

    (options, args) = parser.parse_args(sys.argv)



    ##Change DIR for relative paths:
    base_dir = os.path.split(os.path.abspath(__file__))[0]
    os.chdir(base_dir) #Allow relative running of Rosetta.  Much easier.


    native_db = base_dir+"/"+options.bm_set_db_dir+"/"+"natives."+options.bm_set+"."+ "all"+"."+options.scorefunction+".db3"
    lambda_db = base_dir+"/"+options.bm_set_db_dir+"/"+"natives."+options.bm_set+"."+ "lambda"+"."+options.scorefunction+".db3"
    kappa_db = base_dir+"/"+options.bm_set_db_dir+"/"+"natives."+options.bm_set+"."+ "kappa"+"."+options.scorefunction+".db3"

    analyzer = AnalyzeBenchmarks(native_db, lambda_db, kappa_db)
    analyzer.set_check_nstruct(options.check_nstruct)
    analyzer.set_nstruct_tolerance(options.nstruct_tolerance)
    analyzer.set_delete_dbs(options.delete_dbs)
    analyzer.set_rosetta_extension(options.rosetta_extension)


    BM_info = BenchmarkInfo()
    BM_info.set_data(options.full_name, options.final_name, options.scorefunction, options.decoy_dir, options.exp_batch, options.log_dir, None, options.expected_nstruct)
    analyzer.set_benchmarks([BM_info])




    analyzer.analyze_benchmarks(options.get_ensemble_data, options.skip_antibody_features, options.skip_cluster_features)

    print "Complete"


########################################################################################################################
### Core Classes
########################################################################################################################

class AnalyzeBenchmarks:
    def __init__(self, native_db, lambda_db, kappa_db, benchmarks = [], check_nstruct = False, nstruct_tolerance = 5):
        """
        full_name: $exp.$pdbs.$type.$benchmark
        final_name: docked_relax_top
        """

        self.base_dir = os.path.split(os.path.abspath(__file__))[0]
        self.recovery_base = self.base_dir+"/features/recovery"
        self.ab_db = os.getenv("ROSETTA3_DB")+"/sampling/antibodies/antibody_database_rosetta_design.db"
        self.native_db = native_db
        self.lambda_db = lambda_db
        self.kappa_db = kappa_db
        self.benchmarks = benchmarks

        if not os.path.exists(native_db):
            sys.exit("Native DB path does not exist: "+native_db)

        #Figure out date
        self.date = tools.get_today()
        self.check_nstruct = check_nstruct
        self.nstruct_tolerance = nstruct_tolerance

        self.delete_dbs = True
        self.extension = "macosclangrelease"

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

        if benchmark.has_log_path() and not os.path.exists(benchmark.rel_log_path):
            print "Relative log dir given, but does not exist!"
            print benchmark.rel_log_path
            return False

        if not self.check_nstruct:
            return True

        log_file = "benchmark_nstruct_log.txt"
        today = datetime.date.today()
        today_str = today.strftime("%Y/%m/%d")
        all_found_bool = False

        files = glob(benchmark.decoy_dir+"/"+benchmark.full_name+"*.pdb.gz")
        found_nstruct = len(files)



        if not os.path.exists(log_file):
            FILE = open(log_file, 'w')
            FILE.write("#date name all_clear expected found\n")
        else:
            FILE = open(log_file, 'a')


        if found_nstruct < benchmark.expected_nstruct - self.nstruct_tolerance:
            print "Expected "+repr(benchmark.expected_nstruct)+" decoys but only "+repr(len(files))+" ran."
            FILE.write(today_str + "\t"+benchmark.full_name+"\tFalse\t"+repr(benchmark.expected_nstruct)+"\t"+repr(len(files))+"\n")
            all_found_bool = False
        else:
            FILE.write(today_str + "\t"+benchmark.full_name+"\tTrue\t"+repr(benchmark.expected_nstruct)+"\t"+repr(len(files))+"\n")
            all_found_bool = True

        FILE.close()

        return all_found_bool

    def add_benchmark(self, benchmark):
        if isinstance(benchmark, BenchmarkInfo): pass
        self.benchmarks.append(benchmark)

    def set_benchmarks(self, benchmarks):
        self.benchmarks = benchmarks

    def set_check_nstruct(self, check_nstruct):
        self.check_nstruct = check_nstruct

    def set_nstruct_tolerance(self, nstruct_tolerance):
        self.nstruct_tolerance = nstruct_tolerance

    def set_delete_dbs(self, delete_dbs):
        """
        Set a boolean to delete the present dbs.
        If True, we delete and rerun.
        If False, we keep and skip.
        """
        self.delete_dbs = delete_dbs

    def set_rosetta_extension(self, rosetta_extension):
        self.extension = rosetta_extension

    def analyze_benchmarks(self, get_ensemble_data = False, skip_antibody_features = False, skip_cluster_features = False):
        """
        Runs the FUll analysis
        """
        benchmarks = self._setup_benchmarks()

        for benchmark in benchmarks:
            self._run_analysis(benchmark, get_ensemble_data, skip_antibody_features, skip_cluster_features)

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

            elif self.delete_dbs:
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

                    new_benchmark = copy.deep_copy(benchmark)
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

        r_cmd = "compare_sample_sources.R --config "+json_file
        print r_cmd

        new_recovery_path = self.recovery_base+"/"+benchmark.final_name
        if os.path.exists(new_recovery_path):
            os.system("rm -r "+new_recovery_path)

        os.system(r_cmd)
        os.system("mkdir "+new_recovery_path)
        #os.system(r_cmd)

        #Move cdr cluster features since I still can't seem to fix giving an output dir.
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
            recovery_db_name = benchmark.full_name+"top.recoveries."+benchmark.scorefunction+".db"

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
    """

    Use:
     Creates AntibodyFeatures and ClusterFeatures databases.
     Runs recovery features scripts.
     Adds recovery and risk ratio data to these databases.

    DO NOT ATTEMPT TO MAKE A GUI FOR THIS!


    Example cmd-lines:


                                 SINGLE: RUN ON A SINGLE BENCHMARK

    ./AnalyzeBenchmark.py --mode single --decoy_dir decoys/baseline_hb_idealize_bug_fixes --full_name baseline_hb_idealize_bug_fixes.natives.pack.normal_graft_des_ --final_name pack.normal


    CHECK ALL DECOYS ARE PRESENT

    ./AnalyzeBenchmark.py --mode single --decoy_dir decoys/baseline_hb_idealize_bug_fixes  --full_name baseline_hb_idealize_bug_fixes.natives.pack.normal_graft_des_ --final_name pack.normal --expected_nstruct 200 --check_nstruct


    SPECIFY NATIVE DATABASES

    ./AnalyzeBenchmark.py --mode single --decoy_dir decoys/baseline_hb_idealize_bug_fixes --full_name baseline_hb_idealize_bug_fixes.natives.pack.normal_graft_des_ --final_name pack.normal --native_db xx --kappa_db xx --lambda_db xx


    SPECIFY SCOREFUNCTION USED FOR FEATURES

    ./AnalyzeBenchmark.py --mode single --decoy_dir decoys/baseline_hb_idealize_bug_fixes --full_name baseline_hb_idealize_bug_fixes.natives.pack.normal_graft_des_ --final_name pack.normal --scorefunction orbitals_talaris2013


    """


    main()