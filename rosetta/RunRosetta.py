import json
import os
import sys
import re
import argparse
from RAbD_BM import tools as bm_tools

from rosetta.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral

def run_on_qsub(cmd, queue_dir, name, nodes, ppn, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1

    qsub_cmd = " qsub " +extra_opts+" -V -N "+name + " -d "+queue_dir
    if ppn:

        qsub_cmd = qsub_cmd + " -l nodes="+str(nodes)+":ppn="+str(ppn)
    else:
        qsub_cmd = qsub_cmd + " -l nodes="+str(nodes)

    qsub_cmd = qsub_cmd +" "+script_path
    
    if print_only:
        print(qsub_cmd)
    else:
        os.system(qsub_cmd)


def run_on_slurm(cmd, queue_dir, name, nodes, ppn, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1

    slurm_cmd = "sbatch --nodes="+str(nodes)+" "+extra_opts+" --job-name="+name + " "
    slurm_cmd = slurm_cmd + " -o "+queue_dir+"/"+name+"_%j.out"
    slurm_cmd = slurm_cmd + " -e "+queue_dir+"/"+name+"_%j.err"

    if ppn:
        slurm_cmd = slurm_cmd+" --ntasks="+str(ppn)

    slurm_cmd = slurm_cmd +" "+script_path

    if print_only:
        print(slurm_cmd)
    else:
        os.system(slurm_cmd)

def write_queue_file(cmd, queue_dir, name):
    OUTFILE = open(queue_dir+"/"+name+".sh", 'w')
    OUTFILE.write(cmd)
    OUTFILE.close()
    return queue_dir+"/"+name+".sh"


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

    def _add_args(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--np",
                               default = 101)

        self.parser.add_argument("--nodes",
                                default = 101)

        self.parser.add_argument("--ppn")

        self.parser.add_argument("--nstruct")

        self.parser.add_argument("--platform",
                               default="linux")

        self.parser.add_argument("--compiler", "-c",
                               default = "gcc")


        self.parser.add_argument("--job_manager",
                               default="slurm",
                               help="Job Manager to launch job.",
                               choices = ["slurm","qsub","local"] )

        self.parser.add_argument("--job_manager_opts",
                                 help = "Extra options for the job manager, such as queue or processor requests",
                                 default = "")

        self.parser.add_argument("--machine_file",
                                 help = "Optional machine file for passing to MPI")

        self.parser.add_argument("--print_only",
                                 help = "Do not actually run anything.  Just print out ouput",
                                 default = False,
                                 action = "store_true")


        self.parser.add_argument("--json_base",
                               default = "json_cluster/RAbD_vax.json",
                               help = "JSON file for setting up base paths/etc. for the cluster.")

        self.parser.add_argument("--root",
                                 help = "Override any root directory set in json_base")

        self.parser.add_argument("--job_name",
                                  help = "Override any job name set in json_base")


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
            if hasattr(self.options, "exp_name") and self.options.exp_name:
                pass
            elif self.extra_options.get_exp():
                self.options.exp_name = self.extra_options.get_exp()
            else:
                self.options.exp_name = "unknown_exp"

        def _set_machine_file():
            if self.options.machine_file:
                pass
            elif self.extra_options.get_machine_file():
                self.options.machine_file = self.extra_options.get_machine_file()

            elif self.base_options.get_machine_file():
                self.options.machine_file = self.base_options.get_machine_file()

        def _set_job_manager_opts():
            if self.options.job_manager_opts:
                pass
            elif self.extra_options.get_job_manager_opts():
                self.options.machine_file = self.extra_options.get_job_manager_opts()

            elif self.base_options.get_job_manager_opts():
                self.options.machine_file = self.base_options.get_job_manager_opts()

        #Resolve options overrides
        _set_nstruct()
        _set_exp()
        _set_machine_file()
        _set_job_manager_opts()

    def get_root(self):

        if self.options.root:
            return self.options.root
        elif self.extra_options.get_root():
            return self.extra_options.get_root()
        elif self.base_options.get_root():
            return self.base_options.get_root()
        else:
            return os.getcwd()


    def get_job_name(self, *args, **kwargs):
        """
        Get the job name.
          job_name -> exp_name.

        Pass extra args to add arguments separated by '_'

        """
        if self.options.job_name:

            s = self.options.job_name
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
            os.mkdir(log_path)
        return log_path


    ### Methods to override for specific Benchmarks ###

    def get_make_log_dir(self, *args, **kwargs):
        """
        Get and make the dir to which the MPI output of each process will go.
        """

        print repr(kwargs)
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

        s = s + " -out:prefix "+self._get_out_prefix(*args, **kwargs)


        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self.get_make_out_path(*args, **kwargs)
        s = s + self.base_options.get_base_rosetta_flag_string()

        #Log Dir:
        s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir(*args, **kwargs)

        #For these benchmarks, there are only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        return s

    def get_out_prefix(self, *args, **kwargs):
        """
        Get the output prefix if overrides or exp name is set.
        """
        if self.options.out_prefix:
            return self.options.out_prefix+"."

        s = self.options.exp_name

        return s

    def get_full_cmd(self, *args, **kwargs):
        cmd_string = self.get_output_string(*args, **kwargs)

        if self.options.machine_file:
            cmd = "cd "+ self.get_root()+" \n"+"mpiexec -np " + self.options.np + " --machine_file "+self.options.machine_file+" "+ cmd_string
        else:
            cmd = "cd "+ self.get_root()+" \n"+"mpiexec -np " + self.options.np + " "+ cmd_string

        return cmd


    def run(self, *args, **kwargs):

        print repr(kwargs)
        print "x"
        log_dir = self.get_make_log_dir(*args, **kwargs)
        outpath = self.get_make_out_path(*args, **kwargs)
        queue_dir = self.get_make_queue_dir(*args, **kwargs)

        print "LogDir: "+log_dir
        print "QueueDir: "+queue_dir
        print "OutPath: "+outpath+"\n\n"


        cmd = self.get_full_cmd()

        if self.options.job_manager == "local" and self.options.print_only:
            print cmd + "\n"
        elif self.options.job_manager == "local":
            print cmd + "\n"
            os.system(cmd)

        elif self.options.job_manager == "qsub":
            run_on_qsub(cmd, queue_dir, self.get_job_name(*args, **kwargs), self.options.nodes, self.options.ppn, self.options.print_only, self.options.job_manager_opts)

        elif self.options.job_manager == "slurm":
            run_on_slurm(cmd, queue_dir, self.get_job_name(*args, **kwargs), self.options.nodes, self.options.ppn, self.options.print_only, self.options.job_manager_opts)