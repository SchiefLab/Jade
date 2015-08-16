import json
import os
import sys
import re
import argparse

class RunRosetta(object):
    def __init__(self, program = None):
        """
        Base class for Running Rosetta through python.
        Mainly used for benchmarking experiments.
        Derive this class, override methods to setup benchmark.

        Common Methods to override for more benchmarking control:
          _add_args()
          _get_make_log_dir()
          _get_make_out_path()
          _get_prefix()
          _get_output_string()

          run()

        """
        self.base_options = None
        self.extra_options = None
        self.program = program


        self._add_args()
        self._parse_args()
        self._setup_base_options()
        self._resolve_options()

    def run(self, *args, **kwargs):
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

        self.parser.add_argument("--override_root",
                                 help = "Override the root directory from JSON files")

        self.parser.add_argument("--override_out_prefix")

        self.parser.add_argument("--override_out_path")

        self.parser.add_argument("--override_job_name")

        self.parser.add_argument("--exp_name",
                                 help = "Define the experiment name.  Usually defined in benchmarking for derived classes.")

        if not self.program:
            self.parser.add_argument("--program",
                                 help = "Define the Rosetta program to use.")


    def _parse_args(self):
        self.options = self.parser.parse_args()
        print repr(self.options)
        if hasattr(self.options, 'program'):
            self.program = self.options.program
        elif not self.program:
            sys.exit("Rosetta Program to run must be specified.")


    def _setup_base_options(self):
        """
        Setup the base JSON file that gives settings on the cluster [and] project.
        """
        if not self.options.json_base:
            sys.exit("No Base Json Given.  This is required for general cluster settings.")
        self.base_options = SetupRosettaOptionsGeneral(self.options.json_base)

    def _set_extra_options(self, extra_options):
        """
        Set extra options (derived SetupRosettaOptionsGeneral or baseclass) for benchmarking runs.
        :param extra_options: SetupRosettaOptionsGeneral
        """
        self.extra_options = extra_options
        if not isinstance(extra_options, SetupRosettaOptionsGeneral):
            sys.exit()

    def _resolve_options(self):
        """
        Resolve options from base, extra, and cmd-line settings.
        """
        #Define Conflict resolutions for base class
        def _set_nstruct():
            if self.options.nstruct:
                pass
            elif self.extra_options.get_nstruct():
                self.options.nstruct = str(self.extra_options.get_nstruct())
            elif self.base_options.get_nstruct():
                self.options.nstruct = str(self.base_options.get_nstruct())

        def _set_exp():
            if self.options.exp_name:
                pass
            elif self.extra_options.get_exp():
                self.options.exp_name = self.extra_options.get_exp()
            else:
                self.options.exp_name = "unknown_exp"



        #Resolve options overrides
        _set_nstruct()
        _set_exp()

    def get_root(self):

        if self.options.override_root:
            return self.options.override_root
        elif self.extra_options.get_root():
            return self.extra_options.get_root()
        elif self.base_options.get_root():
            return self.base_options.get_root()
        else:
            return os.getcwd()


    def get_job_name(self, *args):
        """
        Get the job name.
          override_job_name -> exp_name.

        Pass extra args to add arguments separated by '_'

        """
        if self.options.override_job_name:

            s = self.options.override_job_name
            for a in args:
                s = s+"_"+a
            return s

        elif self.extra_options:

            s = self.options.exp_name
            for a in args:
                s = s+"_"+a
            return s

    def _get_program(self):
        """
        Get the set program
        """
        return self.program +".mpi."+self.options.platform + self.options.compiler+"release"

    def get_make_queue_dir(self, *args, **kwargs):
        """
        Get and make the queue dir where qsub/slurm scripts will go.
        """
        log_path = self.base_options.get_make_log_dir()+"/queue_out"
        if not os.path.exists(log_path):
            return log_path


    ### Methods to override for specific Benchmarks ###

    def get_make_log_dir(self, *args, **kwargs):
        """
        Get and make the dir to which the MPI output of each process will go.
        """
        name = self._get_out_prefix(*args, **kwargs)
        log_path = self.base_options.get_make_log_dir()+"/"+name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path



    def get_make_out_path(self, *args, **kwargs):
        """
        Get and make the dir to which decoys will go.  root/decoys
        """
        s = self.base_options.get_root()+"/decoys"
        if not os.path.exists(s):
            os.mkdir(s)

        s = s + "/"+self._get_out_prefix(*args, **kwargs)
        if not os.path.exists(s):
            os.mkdir(s)
        return s


    def get_output_string(self, *args, **kwargs):
        """
        Get the full output string for MPI
        """
        s = self._get_program()

        s = s + " -out:prefix "+self._get_out_prefix(args, kwargs)


        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self.get_make_out_path(args, kwargs)
        s = s + self.base_options.get_base_rosetta_flag_string()

        #Log Dir:
        s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir(args, kwargs)

        #For these benchmarks, there are only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        return s

    def get_out_prefix(self, *args, **kwargs):
        """
        Get the output prefix if overrides or exp name is set.
        """
        if self.options.override_out_prefix:
            return self.options.override_out_prefix+"."

        s = self.options.exp_name

        return s
