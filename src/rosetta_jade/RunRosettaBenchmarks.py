#!/usr/bin/env python

from collections import defaultdict

import argparse
import os
import sys
import json
import re
from collections import defaultdict

from rosetta_jade.SetupRosettaOptionsBenchmark import SetupRosettaOptionsBenchmark
from rosetta_jade.RunRosetta import RunRosetta
from basic import general

from overrides import overrides


class RunRosettaBenchmarks(RunRosetta):
    def __init__(self, program = None, parser = None):
        """
        Derived class for any set of Benchmarks in Rosetta

        """
        RunRosetta.__init__(self, program = program, parser = parser)


        self._current_settings = defaultdict()
        self._current_settings_ordered_keys = []

        self.key_bm_options = 'bm_options'
        self.key_bm_names = 'bm_names'

        self._setup_base_options()

    def _get_list_of_benchmarks(self):
        """
        Get a list of list of the benchmark settings and an ordered list of the benchmark names.

        dict['bm_options'] = list_of_list
        dict['bm_names'] = list of ordered names

        Benchmark names are those that have 'benchmark' and 'rosetta_option' in the list.
        If you would like other benchmark types, sublcass this class and override this function.

        :rtype: defaultdict
        """
        benchmarks = []

        #List the non-rosetta options first, as they are usually more general.
        benchmarks.extend(sorted(self.extra_options.get_non_rosetta_option_benchmark_names()))
        benchmarks.extend(sorted(self.extra_options.get_benchmark_names(only_rosetta=True)))

        list_of_lists = [ self.extra_options.get_benchmarks_of_key(benchmark) for benchmark in benchmarks]

        benchmark_dict = defaultdict()
        benchmark_dict[self.key_bm_options] = list_of_lists
        benchmark_dict[self.key_bm_names] = benchmarks

        if self.options.separate_job_per_pdb and self._get_pdb_list_fname():
            benchmarks.append('pdb')
            list_of_lists.append(self._get_pdb_list_ids())
            self.options.l = None

        self._current_settings_ordered_keys = benchmarks

        return benchmark_dict


    def run_benchmark(self, benchmark_names, benchmark_options):
        """
        Run a single benchmark with options.

        :param benchmark_names: List of benchmark names
        :param benchmark_options: List of benchmark options
        :return:
        """

        for index, bm_name in enumerate(benchmark_names):
            print bm_name+" "+repr(index)
            self._current_settings[bm_name] = benchmark_options[index]

            if bm_name == 'pdb':
                self.options.s = benchmark_options[bm_name]

        self._write_current_benchmark_file()
        RunRosetta.run(self)

    ###################################################################################################################
    ######################################                                #############################################
    #############################                 Full Overrides                    ###################################
    ######################################                                #############################################
    ###################################################################################################################

    @overrides
    def run(self):
        benchmark_dict = self._get_list_of_benchmarks()

        benchmarks_to_run = general.get_all_combos(benchmark_dict[self.key_bm_options])
        for benchmark_set in benchmarks_to_run:
            self.run_benchmark(benchmark_dict[self.key_bm_names], benchmark_set)


    @overrides
    def _add_args(self, parser = None):
        RunRosetta._add_args(self, parser)


        ############################ RAbD Specific Options ################################

        benchmark_options = self.parser.add_argument_group("Benchmark Options", "Options specific for Benchmarking")


        benchmark_options.add_argument("--json_benchmark",
                               help = "JSON file for setting up specific benchmark")

        benchmark_options.add_argument("--pdblist_")
        benchmark_options.add_argument("-p", "--separate_job_per_pdb",
                                default = False,
                                action = "store_true",
                                help = "Separate each PDB in any PDB list given (to python app) into a separate Job and Directory")

    @overrides
    def _setup_base_options(self):
        RunRosetta._setup_base_options(self)

        if not self.options.json_benchmark:
            sys.exit("No Benchmark Json Given.  This is currently required to run benchmarks.")

        self.extra_options = SetupRosettaOptionsBenchmark(self.options.json_benchmark)
        self._resolve_options()

    @overrides
    def _get_make_log_dir(self):

        name = self._get_out_prefix()
        print name
        log_path = self.base_options._get_make_log_dir() + "/" + name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path

    @overrides
    def _get_output_string(self):

        s = RunRosetta._get_output_string(self)


        rosetta_opts = self.extra_options.get_benchmark_names(only_rosetta=True)
        for opt in rosetta_opts:
            print opt
            s = s + self.extra_options.get_rosetta_option_of_key(opt)+" "+str(self._current_settings[opt])

        return s

    @overrides
    def _get_out_prefix(self):
        """
        Setup the prefix from the current options.
        :return:
        """

        json_dict = self.extra_options.json_dict

        if hasattr(self.options, "out_prefix") and self.options.out_prefix:
            return self.options.out_prefix+"."


        s = []

        for key in self._current_settings_ordered_keys:
            if json_dict[key].has_key( self.extra_options.key_use_for_prefix) and not json_dict[key][self.extra_options.key_use_for_prefix]:
                continue

            if self._current_settings[key]:
                opt = self._current_settings[key]
                if type(opt) == bool:
                    if opt == True:
                        opt = "T"
                    else:
                        opt = "F"
                elif type(opt) != str:
                    opt = str(opt)

                s.append(key+"-"+opt)


        return ".".join(s)

    @overrides
    def _get_make_out_path(self):
        """
        Setup the outdir from the current options.
        :return:
        """
        json_dict = self.extra_options.json_dict

        if not os.path.exists("decoys"):
            os.mkdir("decoys")

        if hasattr(self.options, "out_prefix") and self.options.out_prefix:
            return self.options.out_prefix

        s = []

        for key in self._current_settings_ordered_keys:
            if json_dict[key].has_key( self.extra_options.key_use_for_outdir) and not json_dict[key][self.extra_options.key_use_for_outdir]:
                continue

            if self._current_settings[key]:
                opt = self._current_settings[key]
                if type(opt) == bool:
                    if opt == True:
                        opt = "T"
                    else:
                        opt = "F"
                elif type(opt) != str:
                    opt = str(opt)

                s.append(key+"-"+opt)

        outdir = "decoys/"+".".join(s)
        print outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        self._set_outdir(outdir)
        return outdir

    def _get_pdb_list_ids(self):

        pdb_fnames = []
        INFILE = open(self._get_pdb_list_fname(), 'r')
        for line in INFILE:
            line = line.strip()
            if not line or line[0] == "#": continue
            pdb_fnames.append(line)
        INFILE.close()

        return pdb_fnames

    def _get_pdb_list_fname(self):

        if self.options.s:
            return None

        elif not self.options.l and self.options.separate_job_per_pdb:
            sys.exit("PDB LIST must be given for benchmarks using the -l option of the program.  If you are passing it to "
                     "Rosetta proper, pass it to the program instead.")
        else:
            return self.options.l

    def _write_current_benchmark_file(self):
        """
        Writes a file in the directory of decoys which has all the settings that were used for each run.
        :return:
        """

        OUTFILE = open(self._get_make_out_path()+"/RUN_SETTINGS.txt", 'w')
        for benchmark in sorted(self._current_settings):
            OUTFILE.write(" = ".join([benchmark.upper(), str(self._current_settings[benchmark])]) +"\n")
        OUTFILE.close()




if __name__ == "__main__":
    bm = RunRosettaBenchmarks()
    bm.run()