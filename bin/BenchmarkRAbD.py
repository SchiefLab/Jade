#!/usr/bin/env python

from collections import defaultdict

import argparse
import os
import sys
import json
import re
from collections import defaultdict
from rosetta_general.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral
from rosetta_general.RunRosetta import RunRosetta

class AntibodyDesignBMSetup( SetupRosettaOptionsGeneral ):

    def __init__(self, json_file):
        SetupRosettaOptionsGeneral.__init__(self, json_file)

    def get_l_chains(self):
        return self.json_dict["l_chains"]

    def get_exp(self):
        return self.json_dict["exp"]

    def get_outer_cycle_rounds(self):
        if self.json_dict.has_key("outer_cycle_rounds"):
            return self.json_dict["outer_cycle_rounds"]
        else:
            return None

class BenchmarkRAbD(RunRosetta):
    def __init__(self):
        """
        Derived class for running Rosetta Antibody Design (RAbD) Benchmarks
        """
        RunRosetta.__init__(self, "antibody_designer")
        self._set_outer_cycle_rounds()
        self._set_exp()
        self.current_mintype = None
        self.current_l_chain = None

    def _add_args(self):
        RunRosetta._add_args(self)


        ############################ RAbD Specific Options ################################

        self.parser.add_argument("--outer_cycle_rounds")

        self.parser.add_argument("--json_rabd",
                               help = "JSON file for setting up specific benchmark")

        self.parser.add_argument("--mintypes",
                               default = "pack,min,relax")

        self.parser.add_argument("--paper_ab_db",
                               default = True)

        self.parser.add_argument("--with_antigen",
                               default = True)

        self.parser.add_argument("--dataset",
                               default = "new_20")

        self.parser.add_argument("--dock",
                               default = False,
                               action = "store_true")

    def _setup_base_options(self):
        RunRosetta._setup_base_options(self)

        if not self.options.json_rabd:
            sys.exit("No RAbD Json Given.  This is currently required to run benchmarks.")
        self._set_extra_options(AntibodyDesignBMSetup(self.options.json_rabd))

    def _set_outer_cycle_rounds(self):

        if self.options.outer_cycle_rounds:
            pass
        elif self.extra_options.get_outer_cycle_rounds():
            self.options.outer_cycle_rounds = str(self.extra_options.get_outer_cycle_rounds())

    def _set_exp(self):
        if hasattr(self.options, "exp_name") and self.options.exp_name:
            pass
        elif self.extra_options.get_exp():
            self.options.exp_name = self.extra_options.get_exp()
        else:
            self.options.exp_name = "unknown_exp"

    ##Full Overrides###

    def get_make_log_dir(self):

        name = self.get_out_prefix()+"_"+self.current_l_chain
        print name
        log_path = self.base_options.get_make_log_dir()+"/"+name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path

    def get_make_out_path(self):
        s = self.base_options.get_root()+"/decoys"
        if not os.path.exists(s):
            os.mkdir(s)

        s = s + "/"+self.get_out_prefix()
        if not os.path.exists(s):
            os.mkdir(s)
        return s

    def get_output_string(self):
        s = self._get_program()

        s = s + " -out:prefix "+self.get_out_prefix()+"."+" -out:suffix _"+self.current_l_chain+" -antibody:light_chain "+self.current_l_chain


        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self.get_make_out_path()
        s = s + self.base_options.get_base_rosetta_flag_string()

        #Decoys
        s = s + " -l "+self.options.dataset+"."+self.current_l_chain+".PDBLIST.txt"

        #Log Dir:
        if not self.options.one_file_mpi:
            s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir()+"/rosetta_mpi_run"

        #Graft Rounds
        s = s + " -outer_cycle_rounds " + str(self.options.outer_cycle_rounds)

        #Instructions
        s = s + " -cdr_instructions " + self.current_mintype+".instructions.txt"

        #For these benchmarks, there is only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        return s

    def get_out_prefix(self):

        if hasattr(self.options, "out_prefix") and self.options.out_prefix:
            return self.options.out_prefix+"."

        s = ""
        if self.options.with_antigen:
            s = "with_antigen"
        else:
            s = "without_antigen"

        if self.options.dock:
            s = s +".dock"
        else:
            s = s +".no_dock"

        s = s + "."+self.options.exp_name+"."+self.options.dataset


        if self.options.paper_ab_db:
            s = s+"-paper_db"
        else:
            s = s+"-newest_db"

        s = s +"."+self.current_mintype+"."+self.options.outer_cycle_rounds

        return s

    ###Override Run to setup for each mintype and l chain given###
    def run_bm(self):
        mintypesSP = self.options.mintypes.split(",")
        l_chains =  self.extra_options.get_l_chains()

        for mintype in mintypesSP:
            self.current_mintype = mintype
            for l_chain in l_chains:
                self.current_l_chain = l_chain
                self.run()



                """
                cmd_string = ""

                log_dir = self.get_make_log_dir(mintype, l_chain)
                queue_dir = self.get_make_queue_dir()
                outpath = self.get_make_out_path(mintype)

                print "LogDir: "+log_dir
                print "QueueDir: "+queue_dir
                print "OutPath: "+outpath+"\n\n"

                cmd_string = self._get_output_string(mintype, l_chain)
                print cmd_string + "\n"

                if self.options.job_manager == "print_only":
                    return

                if self.options.job_manager == "local":
                    os.chdir(self.base_options.get_root())
                    new_cmd = "cd "+ self.get_root()+" \n"+"mpiexec -np " + self.options.np + " "+ cmd_string
                    print new_cmd
                    #os.system(new_cmd)
                """

if __name__ == "__main__":
    bm = BenchmarkRAbD()
    bm.run_bm()