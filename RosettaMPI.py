import json
import os
import sys
import re
import argparse



class SetupRosettaOptionsGeneral(object):
    """
    Class for setting up more general Rosetta options for benchmarking and repeatable runs on different clusters.  Useful for benchmarking.
    """
    def __init__(self, cluster_json_file):

        FILE = open(cluster_json_file, 'r')
        self.json_dict = json.load(FILE)
        FILE.close()
        self._setup_json_options()

    def get_nstruct(self):
        if self.json_dict.has_key("nstruct"):
            return self.json_dict["nstruct"]
        else:
            return None

    def get_make_log_dir(self):
        log_dir = self.get_root()+"/"+self.json_dict["logs"]
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        return log_dir

    def get_root(self):
        return self._get_root()

    def get_indirs(self):
        return self._get_indirs()

    def get_base_rosetta_flag_string(self, indir_root = None):
        """
        Get the full flag string for output.  Optionally give indir_root for subclasses that require setting of different
         directories, but with same root as given in the cluster file.  Used primarily for complicated benchmarks.

        """
        s = " "

        if len(self._get_flags_files()) > 0:
            for flags_file in self._get_flags_files():
                s = s + " @"+flags_file

        if len(self._get_indirs()) > 0:
            inpaths = " "
            for indir in self._get_indirs():
                if indir_root:
                    inpaths = inpaths+" "+indir_root+"/"+indir
                else:
                    inpaths = inpaths+" "+indir

            s = s+" -in:path "+inpaths

        if len(self._get_set_flags()) > 0:
            for f in self._get_set_flags():
                if not re.search('-', f):
                    f = '-'+f
                s = s+" "+f

        return s

    #################################################
    def _get_root(self):
        """
        Get the root Ouput Dir loaded from JSON. Default is CWD if not set.
        """
        return self.json_dict["root"]

    def _get_flags_files(self):
        """
        Get the list of flags files set in JSON
        """
        return self.json_dict["flags_paths"]

    def _get_set_flags(self):
        return self.json_dict["flags"]

    def _get_indirs(self):
        """
        Get the list of input directories set in JSON
        """
        return self.json_dict["indirs"]

    def _setup_json_options(self):
        if not self.json_dict.has_key("root"):

            self.json_dict["root"] = os.getcwd()
            if not os.path.exists(self._get_root()):
                os.mkdir(self._get_root())

        if not self.json_dict.has_key("flags_files"):

            self.json_dict["flags_files"] = []

        if not self.json_dict.has_key("flags"):
            self.json_dict["flags"] = []

        if not self.json_dict.has_key("indirs"):
            self.json_dict["indirs"] = []

        if not self.json_dict.has_key("flags_paths"):
            self.json_dict["flags_paths"] = []

        if self.json_dict.has_key("indir"):
            for p in self.json_dict["indir"]:
                self.json_dict["indirs"].append(p)

        if self.json_dict.has_key("in_paths"):
            for set in self.json_dict["in_paths"]:
                if not set.has_key("paths"):
                    print "No paths set for in_path.  Skipping"
                    continue

                if set.has_key("root") and set["root"]:
                    root_dir = set["root"]
                    for p in set["paths"]:
                        self.json_dict["indirs"].append(root_dir+"/"+p)
                else:
                    for p in set["paths"]:
                        self.json_dict["indirs"].append(p)

        if self.json_dict.has_key("flags_files"):
            for set in self.json_dict["flags_files"]:
                if not set.has_key("paths"):
                    print "No paths set for flags_files  Skipping"
                    continue

                if set.has_key("root") and set["root"]:
                    root_dir = set["root"]
                    for p in set["paths"]:
                        self.json_dict["flags_paths"].append(root_dir+"/"+p)
                else:
                    for p in set["paths"]:
                        print repr(p)
                        self.json_dict["flags_paths"].append(p)

        ########################################

    def __str__(self):
        return repr(self.json_dict)

    def __repr__(self):
        return repr(self.json_dict)


class RunRosetta(object):
    def __init__(self):

        self.base_options = None
        self.extra_options = None

        self._add_args()
        self._parse_args()
        self._setup_base_options()
        self._resolve_options()

    def run(self, **kwargs):
        pass

    def _add_args(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--np",
                               default = 101)

        self.parser.add_argument("--nstruct")

        self.parser.add_argument("--platform",
                               default="linux")

        self.parser.add_argument("--compiler", "-c",
                               default = "gcc")

        self.parser.add_argument("--job_manager",
                               default="slurm",
                               help="Job Manager to launch job.  Options are: [slurm,qsub,local,print_only]")

        self.parser.add_argument("--json_base",
                               default = "json_cluster/RAbD_vax.json",
                               help = "JSON file for setting up base paths/etc. for the cluster.")

        self.parser.add_argument("--out_prefix_override")

        self.parser.add_argument("--out_path_override")

        self.parser.add_argument("--job_name_override")

        self.parser.add_argument("--exp_name",
                                 help = "Define the experiment name.  Usually defined in benchmarking for derived classes.")

        self.parser.add_argument("--program",
                                 help = "Define the Rosetta program to use.")


    def _parse_args(self):
        self.options = self.parser.parse_args()

    def _setup_base_options(self):
        if not self.options.json_base:
            sys.exit("No Base Json Given.  This is required for general cluster settings.")
        self.base_options = SetupRosettaOptionsGeneral(self.options.json_base)

    def _set_extra_options(self, extra_options):
        self.extra_options = extra_options
        if not isinstance(extra_options, SetupRosettaOptionsGeneral):
            sys.exit()

    def _resolve_options(self):

        #Define Conflict resolutions for base class
        def _set_nstruct():
            if self.options.nstruct:
                pass
            elif self.extra_options.get_nstruct():
                self.options.nstruct = str(self.extra_options.get_nstruct())
            elif self.base_options.get_nstruct():
                self.options.nstruct = str(self.base_options.get_nstruct())

        #Resolve options overrides
        _set_nstruct()

    def get_job_name(self, *args, **kwargs):
        if self.options.job_name_override:
            return self.options.job_name_override
        elif self.extra_options:

            job_name = self.extra_options.get_exp()
            return job_name

    def get_make_log_dir(self, *args, **kwargs):
        name = self._get_out_prefix()+"_"+l_chain
        log_path = self.base_options.get_make_log_dir()+"/"+name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path

    def get_make_queue_dir(self):
        log_path = self.base_options.get_make_log_dir()+"/queue_out"
        if not os.path.exists(log_path):
            return log_path

    def get_make_out_path(self, mintype):
        s = self.base_options.get_root()+"/decoys"
        if not os.path.exists(s):
            os.mkdir(s)

        s = s + "/"+self._get_out_prefix()
        if not os.path.exists(s):
            os.mkdir(s)
        return s


    def get_output_string(self, **kwargs):
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

        #For these benchmarks, there are only a single root directory.
        s = s + self.rabd_options.get_base_rosetta_flag_string(self.base_options.json_dict["in_paths"][0]["root"])

        return s

    def get_out_prefix(self, mintype):
        if self.options.out_prefix_override:
            return self.options.out_prefix_override+"."

        s = ""
        if self.options.with_antigen:
            s = "with_antigen"
        else:
            s = "without_antigen"

        if self.options.dock:
            s = s +".dock"
        else:
            s = s +".no_dock"

        s = s + "."+self.rabd_options.get_exp()+"."+self.options.dataset


        if self.options.paper_ab_db:
            s = s+"-paper_db"
        else:
            s = s+"-newest_db"

        s = s +"."+mintype+"."+self.options.outer_cycle_rounds

        return s

if __name__ == "__main__":
    #For Testing:
    test_file = sys.argv[1]
    setup_mpi = SetupRosettaOptionsGeneral(test_file)
    print repr(setup_mpi)
    print "\n"
    print setup_mpi.get_root()
    print "\n"
    print repr(setup_mpi.get_indirs())
    print "\n"
    print setup_mpi.get_base_rosetta_flag_string()