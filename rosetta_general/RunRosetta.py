#!/usr/bin/env python

import json
import os
import sys
import re
import argparse
from RAbD_BM import tools as bm_tools
from tools.gneral import get_platform

from rosetta_general.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral


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


def run_on_slurm(cmd, queue_dir, name, nodes = False, ntasks = False, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1


    slurm_cmd = "sbatch  "+extra_opts+" --job-name="+name + " "
    slurm_cmd = slurm_cmd + " -o "+queue_dir+"/"+name+"_%j.out"
    slurm_cmd = slurm_cmd + " -e "+queue_dir+"/"+name+"_%j.err"

    if nodes:
        slurm_cmd = slurm_cmd+" --nodes="+str(nodes)
    if ntasks:
        slurm_cmd = slurm_cmd+" --ntasks="+str(ntasks)

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
        if self.options.json_run:
            extra_options = SetupRosettaOptionsGeneral(self.options.json_run)
            self._set_extra_options(extra_options)
        else:
            self._set_extra_options(self.base_options)

        self._resolve_options()

    def _add_args(self):
        self.parser = argparse.ArgumentParser("This program runs Rosetta MPI locally on a cluster using slurm or qsub.")

        if not self.program:
            self.parser.add_argument("--program",
                                 help = "Define the Rosetta program to use if not set in json_run")

        self.parser.add_argument("--np",
                               default = 101)

        self.parser.add_argument("--nodes")

        self.parser.add_argument("--ppn",
                                 help = "Processors per node for qsub.  NTasks is np for slurm")


        self.parser.add_argument("--nstruct")

        self.parser.add_argument("-s",
                                 help = "PDB Files if not set.")

        self.parser.add_argument("-l",
                                 help = "Path to pdb files")

        self.parser.add_argument("--outdir", "-o",
                                 default = "decoys",
                                 help = "Outpath.  "
                                        "Default = 'pwd/decoys' ")

        self.parser.add_argument("--compiler", "-c",
                                 default = "gcc",
                                 help = "Set the compiler used.  Will set clang automatically for macos."
                                        "Default = 'gcc' ",
                                 choices = ["gcc", "clang"])


        self.parser.add_argument("--job_manager",
                               default="slurm",
                               help="Job Manager to launch job. "
                                    "Default = 'slurm ' ",
                               choices = ["slurm","qsub","local"] )


        self.parser.add_argument("--job_manager_opts",
                                 help = "Extra options for the job manager, such as queue or processor requests",
                                 default = "")

        self.parser.add_argument("--machine_file",
                                 help = "Optional machine file for passing to MPI")

        self.parser.add_argument("--print_only",
                                 help = "Do not actually run anything.  Just print out run setup",
                                 default = False,
                                 action = "store_true")


        self.parser.add_argument("--json_base",
                               default = os.path.dirname(os.path.abspath(__file__))+"/jsons/common_flags.json",
                               help = "JSON file for setting up base paths/etc. for the cluster."
                                      "Default = 'file_dir/jsons/common_flags.json' ")

        self.parser.add_argument("--json_run",
                                help = "JSON file for specific Rosetta run.")

        self.parser.add_argument("--root",
                                 help = "Override any root directory set in json_base. If none is set, will use cwd")

        self.parser.add_argument("--job_name",
                                 default = "rosetta_run",
                                help = "Override any job name set in json_base"
                                       "Default = 'rosetta_run'",)


        self.parser.add_argument("--extra_options", "-e",
                                 nargs = '*',
                                 help = "Extra Rosetta options.  "
                                        "Specify like: cdr_instructions=my_file other_option=setting"
                                        "Note NO - charactor. "
                                        "Booleans do not need an = sign.")



    def _parse_args(self):
        self.options = self.parser.parse_args()
        #print repr(self.options)
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

    def _get_extra_rosetta_options_string(self):
        if not self.options.extra_options:
            return ""
        else:
            opts = []
            for o in self.options.extra_options:
                oSP = o.split('=')

                if len(oSP) == 2:
                    opts.append('-'+oSP[0]+" "+oSP[1])
                elif len(oSP) == 1:
                    #Boolean options
                    if oSP[0][0] == '@':
                        opts.append(oSP[0])
                    else:
                        opts.append('-'+oSP[0])
                else:
                    print "Rosetta option too long: "+repr(oSP)
                    continue

            return " ".join(opts)

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

        def _setup_prgram():
            if self.options.program:
                self.program = self.options.program
            elif self.extra_options.get_program():
                self.program = self.extra_options.get_program()
            elif self.base_options.get_prgram():
                self.program = self.base_options.get_program()

        #Resolve options overrides
        _set_nstruct()
        #_set_exp()
        _set_machine_file()
        _set_job_manager_opts()
        _setup_prgram()

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
        else:
            return self.options.exp_name

    def _get_program(self):
        """
        Get the set program
        """

        if get_platform() == 'macos':
            self.options.compiler = 'clang'

        return self.program +".mpi."+get_platform() + self.options.compiler+"release"

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

        #name = self.get_out_prefix(*args, **kwargs)
        log_path = self.base_options.get_make_log_dir()+"/"+self.options.job_name
        if not os.path.exists(log_path):
            os.mkdir(log_path)
        return log_path



    def get_make_out_path(self, *args, **kwargs):
        """
        Get and make the dir to which decoys will go.  root/decoys
        """
        s = self.base_options.get_root()+"/"+self.options.outdir
        if not os.path.exists(s):
            os.mkdir(s)

        #s = s + "/"+self.get_out_prefix(*args, **kwargs)
        if not os.path.exists(s):
            os.mkdir(s)
        return s


    def get_output_string(self, *args, **kwargs):
        """
        Get the full output string for MPI
        """
        s = self._get_program()

        #if self.get_out_prefix(*args, **kwargs):
        #    s = s + " -out:prefix "+self.get_out_prefix(*args, **kwargs)


        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self.get_make_out_path(*args, **kwargs)
        s = s + self.base_options.get_base_rosetta_flag_string()

        #Log Dir:
        s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir(*args, **kwargs)

        #For these benchmarks, there are only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        #Input decoys
        if self.options.l:

            #This works for me - makes it easier to make a list in the dir and run it.
            # Change this to be smarter or add an option if this is a problem.

            p = os.path.dirname(self.options.l)
            s = s+ " -in:path "+p

            s = s+ " -l "+self.options.l

        if self.options.s:
            s = s+ " -s "+self.options.s

        #Extra Rosetta options:
        s = s + " "+self._get_extra_rosetta_options_string()
        return s

    '''
    def get_out_prefix(self, *args, **kwargs):
        """
        Get the output prefix if overrides or exp name is set.
        """
        if self.options.out_prefix:
            return self.options.out_prefix+"."

        else:
            return None

        #s = self.options.exp_name

        #return s
    '''

    def get_full_cmd(self, *args, **kwargs):
        cmd_string = self.get_output_string(*args, **kwargs)

        if self.options.job_manager == "slurm":
            cmd = "cd "+ self.get_root()+" \n"+"mpiexecs "
        else:
            cmd = "cd "+ self.get_root()+" \n"+"mpiexec -np " + self.options.np

        if self.options.machine_file:
            cmd = cmd + " --machine_file "+self.options.machine_file+" "+ cmd_string
        else:
            cmd = cmd + " "+ cmd_string



        return cmd


    def run(self, *args, **kwargs):

        log_dir = self.get_make_log_dir(*args, **kwargs)
        outpath = self.get_make_out_path(*args, **kwargs)
        queue_dir = self.get_make_queue_dir(*args, **kwargs)

        print "\nLogDir: "+log_dir
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
            run_on_slurm(cmd, queue_dir, self.get_job_name(*args, **kwargs), self.options.nodes, self.options.np, self.options.print_only, self.options.job_manager_opts)


if __name__ == "__main__":
    run_rosetta = RunRosetta()
    run_rosetta.run()