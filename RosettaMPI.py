import json
import os
import sys
import re
from optparse import OptionParser, IndentedHelpFormatter



class SetupRosettaOptionsGeneral:
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


class RunRosettaMPI:
    def __init__(self, setup_subclass = None):


        self._parse_args()

    def _parse_args(self):
        parser = OptionParser()

        parser.add_option("--np")
        parser.add_option("--nstruct")
        parser.add_option("--platform",
                          default="linux")

        parser.add_option("--compliler", "-c",
                          default = "gcc")

        parser.add_option("--job_manager",
                          default="slurm")

        parser.add_option("--base_json")

        parser.add_option("--app")

        options, args = parser.parse_args(sys.argv)

    def get_full_app_name(self):
        return self.app + ".mpi."+self.platform + self.compiler+"release"


    def run(self):
        pass
        #Create MPI Job Setup.

        #Create temp script directory in root_outdir

        #Add -separate_outputs


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