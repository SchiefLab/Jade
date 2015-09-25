#!/usr/bin/env python
# Jared Adolf-Bryfogle
# Classes for running Rosetta on a cluster using pre-defined sets of options stored in Json file.

# Used by itself or subclassed for benchmarking.
# Can also use in other code by passing the parser to RunRosetta.

import os
import sys
import re
import argparse
from tools.general import get_platform
from tools.path import *

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
        print(cmd+"\n\n")
        print(qsub_cmd)
    else:
        print(cmd+"\n\n")
        print(qsub_cmd)
        os.system(qsub_cmd)


def run_on_slurm(cmd, queue_dir, name, nodes = False, ntasks = False, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1


    slurm_cmd = "sbatch  "+extra_opts+" --job-name="+name + " "
    slurm_cmd = slurm_cmd + " -o "+queue_dir+"/"+name+"_%j.out"
    slurm_cmd = slurm_cmd + " -e "+queue_dir+"/"+name+"_%j.err"

    if nodes and not re.search('--nodes', extra_opts):
        slurm_cmd = slurm_cmd+" --nodes="+str(nodes)
    if ntasks and not re.search('--ntasks', extra_opts):
        slurm_cmd = slurm_cmd+" --ntasks="+str(ntasks)

    slurm_cmd = slurm_cmd +" "+script_path

    if print_only:
        print cmd + "\n\n"
        print(slurm_cmd)
    else:
        print cmd + "\n\n"
        print(slurm_cmd)
        os.system(slurm_cmd)

def write_queue_file(cmd, queue_dir, name):

    cmd = "#!/bin/bash\n\n\n" +cmd
    OUTFILE = open(queue_dir+"/"+name+".sh", 'w')
    OUTFILE.write(cmd)
    OUTFILE.close()
    return queue_dir+"/"+name+".sh"


class RunRosetta(object):
    def __init__(self, program = None, parser = None):
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


        self._add_args(parser)
        self._parse_args()
        self._setup_base_options()

        if self.options.json_run:
            extra_options = SetupRosettaOptionsGeneral(self.options.json_run)
            self._set_extra_options(extra_options)
        else:
            self._set_extra_options(self.base_options)

        self._resolve_options()

    def _add_args(self, parser = None):
        """
        Add Arguments to an Argument Parser or create a new one.
        """
        if not parser:
            self.parser = argparse.ArgumentParser("This program runs Rosetta MPI locally or on a cluster using slurm or qsub.  "
                                              "Relative paths are accepted.")
        else:
            self.parser = parser



        job_setup = self.parser.add_argument_group("Job Setup" )
        job_setup.add_argument("--job_manager",
                               default="slurm",
                               help="Job Manager to launch job. "
                                    "Default = 'slurm ' ",
                               choices = ["slurm","qsub","local"] )


        job_setup.add_argument("--job_manager_opts",
                                 help = "Extra options for the job manager, such as queue or processor requests"
                                        "Remove double dashes. Specify like: -p imperial exclusive.",
                                 default = [],
                                 nargs = "*")


        job_setup.add_argument("--np",
                               default = 101,
                                 help = "Number of processors to use for MPI.  "
                                        "Default = 101")

        job_setup.add_argument("--nodes",
                                 help = "Number of nodes to ask for.  Optional. ")

        job_setup.add_argument("--ppn",
                                 help = "Processors per node for qsub.  NTasks is np for slurm")



        job_setup.add_argument("--nstruct",
                                 default = 1)

        job_setup.add_argument("--compiler", "-c",
                                 default = "gcc",
                                 help = "Set the compiler used.  Will set clang automatically for macos. "
                                        "Default = 'gcc' ",
                                 choices = ["gcc", "clang"])

        job_setup.add_argument("--machine_file",
                                 help = "Optional machine file for passing to MPI")

        job_setup.add_argument("--job_name",
                                default = "rosetta_run",
                                help = "Set the job name used for mpi_tracer_to_file dir and queue.  "
                                       "Default = 'rosetta_run'.  "
                                       "(Benchmarking: Override any set in json_base.)",)


        protocol_setup = self.parser.add_argument_group("Protocol Setup")

        if not self.program:
            protocol_setup.add_argument("--program",
                                 help = "Define the Rosetta program to use if not set in json_run")

        protocol_setup.add_argument("-s",
                                 help = "Path to a pdb file")

        protocol_setup.add_argument("-l",
                                 help = "Path to a list of pdb files")

        protocol_setup.add_argument("--outdir", "-o",
                                 default = "decoys",
                                 help = "Outpath.  "
                                        "Default = 'pwd/decoys' ")


        protocol_setup.add_argument("--json_base",
                               default = os.path.dirname(os.path.abspath(__file__))+"/jsons/common_flags.json",
                               help = "JSON file for setting up base paths/etc. for the cluster."
                                      "Default = 'file_dir/jsons/common_flags.json' ")

        protocol_setup.add_argument("--json_run",
                                help = "JSON file for specific Rosetta run.  Not required.")

        protocol_setup.add_argument("--root",
                                 help = "Set the root directory.  "
                                        "Default = pwd.  "
                                        "(Benchmarking: Override any set in json_base.)")

        protocol_setup.add_argument("--extra_options", "-e",
                                 nargs = '*',
                                 help = "Extra Rosetta options.  "
                                        "Specify like: cdr_instructions=my_file other_option=setting.  "
                                        "Note NO - charactor. "
                                        "Booleans do not need an = sign.")

        protocol_setup.add_argument("--one_file_mpi",
                                 help = "Don't setup mpi_tracer_to_file. ",
                                 default = False,
                                 action = "store_true")

        protocol_setup.add_argument("--print_only",
                                 help = "Do not actually run anything.  Just print setup for review.",
                                 default = False,
                                 action = "store_true")

        db_group = self.parser.add_argument_group("Relational Databases", "Options for Rosetta Database input and output.  Use for features or for inputting and output structures as databases")

        db_group.add_argument("--db_mode",
                            help = "Set the mode for Rosetta to use if using a database.  "
                                      "Features will be output to a database.  "
                                      "If not sqlite3, must build Rosetta with extras.  "
                                      "If any post-processing is required, such as combining sqlite3 dbs, will do this.  "
                                      "Default DB mode for features is sqlite3.  ",
                            choices = ["sqlite3", "mysql", "postgres"])

        db_group.add_argument("--db_name",
                              help = "In or Out database name",
                              default = "features.db")

        db_group.add_argument("--db_batch",
                              help = "Batch of structures.",
                              default ="feat")

        db_group.add_argument("--db_in",
                              help = "Use an input database",
                              default = False,
                              action = "store_true")

        db_group.add_argument("--db_out",
                              help = "Use an output database",
                              default = False,
                              action = "store_true")

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

    def set_json_run(self, json_run):
        self._set_extra_options(SetupRosettaOptionsGeneral( json_run))
    def _set_extra_options(self, extra_options):
        """
        Set extra options (derived SetupRosettaOptionsGeneral or baseclass) for benchmarking runs.
        :param extra_options: SetupRosettaOptionsGeneral
        """
        self.extra_options = extra_options
        if not isinstance(extra_options, SetupRosettaOptionsGeneral):
            sys.exit()
        self._resolve_options()

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

        def _setup_program():
            if hasattr(self.options, "program") and self.options.program:
                self.program = self.options.program
            elif self.extra_options.get_program():
                self.program = self.extra_options.get_program()
            elif self.base_options.get_program():
                self.program = self.base_options.get_program()

        def _set_db_mode():
            if self.options.db_mode:
                pass
            elif self.extra_options.get_db_mode():
                self.options.db_mode = self.extra_options.get_db_mode()
            elif self.base_options.get_db_mode():
                self.options.db_mode = self.base_options.get_db_mode()

        def _clean_up_db_name():
            if not re.search(".db", self.options.db_name):
                self.options.db_name = self.options.db_name+".db"

        #Resolve options overrides
        _set_nstruct()
        #_set_exp()
        _set_db_mode()
        _clean_up_db_name()
        _set_machine_file()
        _set_job_manager_opts()
        _setup_program()


    def get_root(self):

        if self.options.root:
            return self.options.root
        elif self.extra_options.get_root():
            return self.extra_options.get_root()
        elif self.base_options.get_root():
            return self.base_options.get_root()
        else:
            return os.getcwd()

    def get_job_manager_opts(self):
        opts = []
        for opt in self.options.job_manager_opts:
            if re.search('-', opt):
                opts.append(opt)
            else:
                opts.append("--"+opt)
        return " ".join(opts)

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
        log_path = self.base_options.get_make_log_dir()+"/queue"
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


    def set_outdir(self, outdir):
        self.options.outdir = outdir
        self.get_make_out_path()

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
        if not self.options.one_file_mpi:
            s = s + " -mpi_tracer_to_file "+ self.get_make_log_dir(*args, **kwargs)+"/rosetta_mpi_run"

        #For these benchmarks, there are only a single root directory.
        s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options.get_root())

        #DB Mode
        if self.options.db_in:
            s = s + " -in:use_database"
            if not self.options.db_mode:
                sys.exit("Please select the database mode you wish to use.")

        if self.options.db_out:
            s = s + " -out:use_database"
            if not self.options.db_mode:
                sys.exit("Please select the database mode you wish to use. ")

            if self.options.db_mode == "sqlite3" and (not re.search('separate_db_per_mpi_process', s)) :
                s = s + " -separate_db_per_mpi_process"

        if self.options.db_mode:
            s = s + " -inout:dbms:mode "+self.options.db_mode

        if self.options.db_name:
            s = s + " -inout:dbms:database_name " +self.options.db_name


        if re.search("features", self.base_options.get_xml_script() + self.extra_options.get_xml_script()):

            if self.options.db_name:
                s = s +" -parser:script_vars name="+self.options.db_name
            if self.options.db_batch:
                s = s +" -parser:script_vars batch="+self.options.db_batch


        #Input decoys
        if self.options.l:

            #This works for me - makes it easier to make a list in the dir and run it.
            # Change this to be smarter or add an option if this is a problem.

            p = os.path.dirname(self.options.l)
            if p != "":
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
            cmd = "cd "+ self.get_root()+" \n"+"mpiexec "
        else:
            cmd = "cd "+ self.get_root()+" \n"+"mpiexec -np " + str(self.options.np)

        if self.options.machine_file:
            cmd = cmd + " --machine_file "+self.options.machine_file+" "+ cmd_string
        else:
            cmd = cmd + " "+ cmd_string

        if self.options.db_mode == "sqlite3":
            cmd = cmd + "\n"
            cmd = cmd + "cd "+self.options.outdir+"\n"
            cmd = cmd + "bash "+ get_rosetta_features_root()+"/sample_sources/merge.sh "+self.options.db_name + " "+self.options.db_name+"_*\n"
            cmd = cmd + "cd -"

        if re.search("out:path:all", cmd):
            sys.exit( "Please use --outdir script option instead of out:path:all" )
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
            run_on_qsub(cmd, queue_dir, self.get_job_name(*args, **kwargs), self.options.nodes, self.options.ppn, self.options.print_only, self.get_job_manager_opts())

        elif self.options.job_manager == "slurm":
            run_on_slurm(cmd, queue_dir, self.get_job_name(*args, **kwargs), self.options.nodes, self.options.np, self.options.print_only, self.get_job_manager_opts())


if __name__ == "__main__":
    run_rosetta = RunRosetta()
    run_rosetta.run()