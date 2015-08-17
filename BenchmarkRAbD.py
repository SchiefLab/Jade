from collections import defaultdict

import argparse
import os
import sys
import json
import re
from collections import defaultdict
from rosetta.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral
from rosetta.RunRosetta import RunRosetta

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

        print "Where are my options" + repr(self.options)


        if self.options.outer_cycle_rounds:
            pass
        elif self.extra_options.get_outer_cycle_rounds():
            self.options.outer_cycle_rounds = str(self.extra_options.get_outer_cycle_rounds())


    ##Full Overrides###

    def get_make_log_dir(self, mintype, l_chain):
        name = self._get_out_prefix(mintype)+"_"+l_chain
        log_path = self.base_options.get_make_log_dir()+"/"+name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path

    def get_make_out_path(self, mintype):
        s = self.base_options.get_root()+"/decoys"
        if not os.path.exists(s):
            os.mkdir(s)

        s = s + "/"+self._get_out_prefix(mintype)
        if not os.path.exists(s):
            os.mkdir(s)
        return s

    def _get_output_string(self, mintype, l_chain):
        s = self._get_program()

        s = s + " -out:prefix "+self._get_out_prefix(mintype)+"."+" -out:suffix _"+l_chain+" -antibody:light_chain "+l_chain


        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self.get_make_out_path(mintype)
        s = s + self.base_options.get_base_rosetta_flag_string()

        #Decoys
        s = s + " -l "+self.options.dataset+"."+l_chain+".PDBLIST.txt"

        #Log Dir:
        s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir(mintype, l_chain)

        #Graft Rounds
        s = s + " -outer_cycle_rounds " + str(self.options.outer_cycle_rounds)

        #Instructions
        s = s + " -cdr_instructions " + mintype+".instructions.txt"

        #For these benchmarks, there is only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        return s

    def _get_out_prefix(self, mintype):
        if self.options.override_out_prefix:
            return self.options.override_out_prefix+"."

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

        s = s +"."+mintype+"."+self.options.outer_cycle_rounds

        return s

    ###Override Run to setup for each mintype and l chain given###
    def run(self):
        mintypesSP = self.options.mintypes.split(",")
        l_chains =  self.extra_options.get_l_chains()

        for mintype in mintypesSP:
            for l_chain in l_chains:
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
                    new_cmd = "cd "+ self.base_options.get_root()+" \n"+"mpiexec -np " + self.options.np + " "+ cmd_string
                    print new_cmd
                    #os.system(new_cmd)

if __name__ == "__main__":
    bm = BenchmarkRAbD()
    bm.run()