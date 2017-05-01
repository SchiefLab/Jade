import os
import sys
from jade.rosetta_jade.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral





class SetupRosettaOptionsBenchmark(SetupRosettaOptionsGeneral):
    """
    Class for setting up Rosetta Benchmarks.  See database/rosetta/benchmark_jsons_rabd/nstruct_test.json for an example.

    Basically, a set of benchmarks and rosetta options are given in the JSON.
    Other keys can be specified for specific benchmarks (like the instructions file stuff in the above file.)

    This can be used to use a single JSON file and run RosettaMPI on ALL combinations of benchmarks given.
    """
    def __init__(self, json_file):
        SetupRosettaOptionsGeneral.__init__( self, json_file)

        self.key_benchmarks = "benchmarks"
        self.key_rosetta_option = "rosetta_option"

        self.key_use_for_prefix = "use_for_prefix"
        self.key_use_for_outdir = "use_for_outdir"

    def get_exp(self):
        """
        Get the benchmark name or fail.
        :rtype: str
        """
        if self.json_dict.has_key("exp"):
            return self.json_dict["exp"]
        else:
            sys.exit("Please include the key 'exp' to define overall the name of the benchmarking experiment.")

    def get_benchmark_names(self, only_rosetta = False):
        """
        Get the names of all the benchmarks we will run.

        Each benchmark must have a dictionary that defines 'benchmarks' as a list.
        You may optionally give the rosetta_option.
        Currently, your subclass of RunRosetta will need to code how all this is run.  Hopefully, that will change.

        If only_rosetta is true, will only give the benchmark names that are based on rosetta options.

        For example:

        "outer_cycle_rounds":{
            "rosetta_option":"-outer_cycle_rounds",
            "benchmarks":[ 25, 50, 75, 100]
        },

        :rtype: list
        """

        benchmark_names = []
        for key in self.json_dict:
            if type(self.json_dict[key]) == dict and self.json_dict[key].has_key(self.key_benchmarks):
                if only_rosetta and type(self.json_dict[key]) == dict and self.json_dict[key].has_key(self.key_rosetta_option):
                    benchmark_names.append(key)
                elif not only_rosetta:
                    benchmark_names.append(key)
                else:
                    continue

        return benchmark_names

    def get_non_rosetta_option_benchmark_names(self):
        """
        Similar to get_benchmark_names, but only for options which do not have the tag rosetta_option

        :rtype: list
        """

        benchmark_names = []
        for key in self.json_dict:
            if type(self.json_dict[key]) == dict and self.json_dict[key].has_key(self.key_benchmarks):
                if  not self.json_dict[key].has_key(self.key_rosetta_option):
                    benchmark_names.append(key)
                else:
                    continue

        return benchmark_names


    def get_benchmarks_of_key(self, benchmark_name):
        """
        Get the list of benchmarks for a particular benchmark key.
        :param benchmark_name: str
        :rtype: list
        """
        if not self.json_dict.has_key(benchmark_name):
            sys.exit("Could not find benchmark name in json dict!  "+benchmark_name)
        else:

            #In special circumstances, it may be a list (as for the CDRs), since they can both be a benchmark or all together.

            try:
                return self.json_dict[benchmark_name][self.key_benchmarks]
            except TypeError:
                return self.json_dict[benchmark_name]

    def get_rosetta_option_of_key(self, benchmark_name):
        """
        Get the Rosetta option
        :param benchmark_name:
        :rtype: str
        """

        return self.json_dict[benchmark_name][self.key_rosetta_option]

    def use_benchmark_for_outdir(self, benchmark):
        """
        Should we use the benchmark name for output?

         Specified by the 'use_for_outdir' in JSON.
         If not specified, or benchmark not in list, we assume True!

        :param benchmark: str
        :rtype: bool
        """


        if self.json_dict.has_key(benchmark):
            if self.json_dict[benchmark].has_key(self.key_use_for_outdir):
                if not self.json_dict[benchmark][self.key_use_for_outdir]:
                    return False

        return True

    def use_benchmark_for_prefix(self, benchmark):
        """
        Should we use the benchmark name for prefix?

         Specified by the 'use_for_prefix' in JSON.
         If not specified, or benchmark not in list, we assume True!

        :param benchmark: str
        :rtype: bool
        """


        if self.json_dict.has_key(benchmark):
            if self.json_dict[benchmark].has_key(self.key_use_for_prefix):
                if not self.json_dict[benchmark][self.key_use_for_prefix]:
                    return False

        return True