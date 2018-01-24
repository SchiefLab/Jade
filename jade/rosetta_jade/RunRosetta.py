#!/usr/bin/env python
# Jared Adolf-Bryfogle
# Classes for running Rosetta on a cluster using pre-defined sets of options stored in Json file.

# Used by itself or subclassed for benchmarking.
# Can also use in other code by passing the parser to RunRosetta.

import os
import sys
from collections import defaultdict
import argparse
from jade.basic.general import get_platform
from jade.basic.path import *
from jade.basic.general import fix_input_args

from jade.rosetta_jade.SetupRosettaOptionsGeneral import SetupRosettaOptionsGeneral

#Fixes parser for extra rosetta opts.
fix_input_args()

def run_on_qsub(cmd, queue_dir, name, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1

    qsub_cmd = " qsub " +extra_opts+" -V -N "+name + " -d "+queue_dir

    qsub_cmd = qsub_cmd +" "+script_path

    if print_only:
        print print_full_cmd(cmd, script_path)
        print("\n\n")
        print(qsub_cmd)
    else:
        #qsub_cmd = "which sbatch"
        print_full_cmd(cmd, script_path)

        print "\n\n"

        print(qsub_cmd)
        os.system(qsub_cmd)


def run_on_slurm(cmd, queue_dir, name, nodes = False, ntasks = None, print_only = False, extra_opts = ""):
    script_path = write_queue_file(cmd, queue_dir, name)


    #qsub -q dna -l nodes=10:ppn=11 -V -N $1 -d $qsub_output -v np=101 $benchmarks/$1

    ##Set the walltime for something rediculous.  I hate walltime.
    extra_opts = extra_opts+" "+"--time=0"

    slurm_cmd = "sbatch  "+extra_opts+" --job-name="+name + " "
    slurm_cmd = slurm_cmd + " -o "+queue_dir+"/"+name+"_%j.out"
    slurm_cmd = slurm_cmd + " -e "+queue_dir+"/"+name+"_%j.err"

    if nodes and not re.search('--nodes', extra_opts):
        slurm_cmd = slurm_cmd+" --nodes="+str(nodes)
    if ntasks and not re.search('--ntasks', extra_opts):
        slurm_cmd = slurm_cmd+" --ntasks="+str(ntasks)

    slurm_cmd = slurm_cmd +" "+script_path

    if print_only:
        print "Only Printing!"

        print_full_cmd(cmd, script_path)

        print "\n\n"
        print(slurm_cmd)
    else:
        print "Starting Slurm!!"
        
        #slurm_cmd = "which sbatch"

        print_full_cmd(cmd, script_path)
        os.system("which sbatch")
        os.system(slurm_cmd)

def print_full_cmd(cmd, script_path = None):
    cmdSP = cmd.split(" ")
    print cmdSP[0]+" "+cmdSP[1]
    print "\n"

    flags = get_option_strings(cmd)

    print flags

    if script_path:
        outfilename = script_path.replace(".sh", ".flags")

        F = open(outfilename, 'w')
        F.write(flags)
        F.close()
        print "\nFlags file written to "+outfilename


def get_option_strings(cmd):
    """
    Get the options as a string to be printed or saved to a file.
    :param cmd:
    :rtype: str
    """

    grouped = defaultdict(list)
    cmdSP = cmd.split(" ")
    if len(cmdSP)< 3: return ""
    options = []
    current_option = cmdSP[2]

    for c in cmdSP[3:]:
        if not c: continue
        if c[0] =='-' and c != '-':
            current_option = c
            options.append(current_option)

            continue
        elif c != '-':
            grouped[current_option].append(c)

    final_string = []
    for option in options:
        opt = option+" "+" ".join(grouped[option])
        if opt[0] =='-' and opt != '-':
            final_string.append(opt)

    full_s = "\n".join(final_string)


    return "\n".join([s for s in full_s.split("\n") if (s[0] == '-' and s  != '-')])


def write_queue_file(cmd, queue_dir, name):

    cmd = "#!"+os.environ['SHELL']+"\n\n\n" +cmd
    OUTFILE = open(queue_dir+"/"+name+".sh", 'w')
    OUTFILE.write(cmd)
    OUTFILE.close()
    return queue_dir+"/"+name+".sh"


class RunRosetta(object):
    def __init__(self, program = None, parser = None, db_mode = False, json_run=None):
        """
        Base class for Running Rosetta through python.
        Mainly used for benchmarking experiments.
        Derive this class, override methods to setup benchmark.

        Common Methods to override for more benchmarking control:
          _add_args()
          _get_make_mpi_tracer_dir()
          _get_make_out_path()
          _get_prefix()
          _get_output_string()

          run()

        """
        self.base_options = None
        self.extra_options = None
        self.program = program
        self.db_mode = db_mode

        self.jsons = [os.path.basename(d) for d in glob.glob(os.path.join(get_rosetta_json_run_path(),"*.json"))]

        self._add_args(parser)
        self._parse_args()
        self._setup_base_options()

        if json_run:
            self.options.json_run = json_run

        if self.options.json_run:
            print "Setting JSON RUN"
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


        common_options = self.parser.add_argument_group("Common Options" )

        common_options.add_argument("-s",
                                 help = "Path to a pdb file")

        common_options.add_argument("-l",
                                 help = "Path to a list of pdb files")

        common_options.add_argument("--np",
                               default = 101,
                                 help = "Number of processors to use for MPI.  "
                                        "Default = 101")

        common_options.add_argument("--nstruct",
                               default = 1,
                               help = "The number of structures/parallel runs.  Can also set this in any JSON file.")


        common_options.add_argument("--job_name",
                                default = "rosetta_run",
                                help = "Set the job name used for mpi_tracer_to_file dir and queue.  "
                                       "Default = 'rosetta_run'.  "
                                       "(Benchmarking: Override any set in json_base.)",)

        common_options.add_argument("--outdir", "-o",
                                 default = "decoys",
                                 help = "Outpath.  "
                                        "Default = 'pwd/decoys' ")

        common_options.add_argument("--json_run",
                                help = "JSON file for specific Rosetta run.  Not required.  Pre-Configured JSONS include: "+repr(self.jsons),
                                )


        common_options.add_argument("--extra_options",
                                 help = "Extra Rosetta options.  "
                                        "Specify in quotes!")

        common_options.add_argument("--script_vars",
                                help = "Any script vars for XML scripts."
                                       "Specify as you would in Rosetta. like: glycosylation=137A,136A",
                                nargs = '*')

        common_options.add_argument("--jd3", help = "Is this app JD3?  Must build with extras=mpi,serialization.",
                                    default = False,
                                    action="store_true")

        if not self.program:
            common_options.add_argument("--program",
                                 help = "Define the Rosetta program to use if not set in json_run")




        debug_options = self.parser.add_argument_group("Testing and Debugging")

        debug_options.add_argument("--print_only",
                                 help = "Do not actually run anything.  Just print setup for review.",
                                 default = False,
                                 action = "store_true")


        debug_options.add_argument("--local_test",
                               default = False,
                               help = "Is this a local test?  Will change nstruct to 1 and run on 2 processors",
                               action = "store_true")

        debug_options.add_argument("--one_file_mpi",
                                 help = "Output all MPI std::out to a single file instead of splitting it. ",
                                 default = False,
                                 action = "store_true")



        special_options = self.parser.add_argument_group("Special Options for controlling execution")
        special_options.add_argument("--job_manager",
                               default="slurm",
                               help="Job Manager to launch job. (Or none if local or local_test)"
                                    "Default = 'slurm ' ",
                               choices = ["slurm","qsub","local","local_test"] )



        special_options.add_argument("--job_manager_opts",
                                 help = "Extra options for the job manager, such as queue or processor requests"
                                        "Remove double dashes. Exclusive is on by default.  Specify like: -p imperial exclusive.",
                                 default = [],
                                 nargs = "*")

        special_options.add_argument("--json_base",
                               default = get_rosetta_json_run_path()+"/common_flags.json",
                               help = "JSON file for setting up base paths/etc. for the cluster."
                                      "Default = 'database/rosetta/jsons/common_flags.json' ")

        special_options.add_argument("--compiler",
                                 default = "gcc",
                                 help = "Set the compiler used.  Will set clang automatically for macos. "
                                        "Default = 'gcc' ",
                                 choices = ["gcc", "clang"])




        special_options.add_argument("--mpiexec",
                               help = "Specify a particular path (or type of) MPI exec. Default is srun (due to vax). If local or local test, will use mpiexex",
                               default = "srun")

        special_options.add_argument("--machine_file",
                                 help = "Optional machine file for passing to MPI")








        if self.db_mode:
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
                              default = "features.db"
                                )

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

        if self.options.local_test:
            self.options.job_manager = "local_test"

    def _setup_base_options(self):
        """
        Setup the base JSON file that gives settings on the cluster [and] project.
        """
        if not self.options.json_base:
            sys.exit("No Base Json Given.  This is required for general cluster settings.")
        self.base_options = SetupRosettaOptionsGeneral(self.options.json_base)

    def _set_json_run(self, json_run):
        print "JSON Run Set: "+json_run
        self.options.json_run = json_run
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

            '''
            opts = []

            skip_next = False
            for i in range(0, len(self.options.extra_options)):

                o = self.options.extra_options[ i ]
                if skip_next:
                    skip_next = False
                    continue

                oSP = o.split('=')

                if len(oSP) == 2:
                    opts.append('-'+oSP[0]+" "+oSP[1])
                elif len(oSP) == 1:
                    #Boolean options
                    if oSP[0][0] == '@' and len(oSP[0]) != 1:
                        opts.append(oSP[0])
                    elif oSP[0][0] == '@':
                        opts.append(oSP[0])

                        opts.append(self.options.extra_options[ i + 1 ])
                        skip_next = True

                    else:
                        opts.append('-'+oSP[0])
                else:
                    print "Rosetta option too long: "+repr(oSP)
                    continue

            return " ".join(opts)
            '''

            return self.options.extra_options

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
            elif self.extra_options._get_job_manager_opts():
                self.options.machine_file = self.extra_options._get_job_manager_opts()

            elif self.base_options._get_job_manager_opts():
                self.options.machine_file = self.base_options._get_job_manager_opts()

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
            if self.options.db_name:
                if not re.search(".db", self.options.db_name):
                    self.options.db_name = self.options.db_name+".db"

        #Resolve options overrides
        _set_nstruct()
        _set_machine_file()
        _set_job_manager_opts()
        _setup_program()

        if self.db_mode:
            _set_db_mode()
            _clean_up_db_name()

    def _get_root(self):

        if self.extra_options._get_root():
            return self.extra_options._get_root()
        elif self.base_options._get_root():
            return self.base_options._get_root()
        else:
            return os.getcwd()

    def _get_job_manager_opts(self):
        opts = []
        for opt in self.options.job_manager_opts:
            if re.search('-', opt):
                opts.append(opt)
            else:
                opts.append("--"+opt)
        return " ".join(opts)

    def _get_job_name(self, *args, **kwargs):
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
            return self.options.exp_name+"."+os.path.basename(self._get_make_out_path())

    def _get_program(self):
        """
        Get the set program
        """

        if get_platform() == 'macos':
            self.options.compiler = 'clang'

        if not self.program and not self.options.json_run:
            sys.exit("Please set a JSON run file.")

        elif not self.program:
            sys.exit("Please specify a program in the JSON run file")

        if self.options.jd3:
            mode=".mpiserialization."
        else:
            mode=".mpi."

        return self.program +mode+get_platform() + self.options.compiler+"release"

    def _get_make_queue_dir(self, *args, **kwargs):
        """
        Get and make the queue dir where qsub/slurm scripts will go.
        """
        log_path = self.base_options._get_make_log_root_dir() + "/queue"
        if not os.path.exists(log_path):
            os.makedirs(log_path)
        return log_path

    def _set_outdir(self, outdir):
        self.options.outdir = outdir


    ### Methods to override for specific Benchmarks ###

    def _get_make_mpi_tracer_dir(self, *args, **kwargs):
        """
        Get and make the dir to which the MPI output of each process will go.
        ONLY for MPI TRACER LOGs
        """

        #name = self.get_out_prefix(*args, **kwargs)
        log_path = os.path.join(self.base_options._get_make_log_root_dir(), "mpi_tracer_logs")
        if not os.path.exists(log_path):
            os.makedirs(log_path)

        log_path = log_path+"/"+self._get_job_name()
        if not os.path.exists(log_path):
            os.makedirs(log_path)
        return log_path

    def _get_out_prefix(self, *args, **kwargs):
        return None

    def _get_make_out_path(self, *args, **kwargs):
        """
        Get and make the dir to which decoys will go.  root/decoys
        """
        s = self.base_options._get_root() + "/" + self.options.outdir
        if not os.path.exists(s):
            os.makedirs(s)
        #s = s + "/"+self.get_out_prefix(*args, **kwargs)
        if not os.path.exists(s):
            os.makedirs(s)
        return s



    def _get_output_string(self, *args, **kwargs):
        """
        Get the full output string for MPI
        """
        s = self._get_program()

        if self._get_out_prefix():
            s = s + " -out:prefix "+self._get_out_prefix()

        if re.search("out:path:all", self._get_extra_rosetta_options_string()):
            #print self._get_extra_rosetta_options_string()
            sys.exit( "Please use --outdir script option instead of out:path:all" )



        #Nstruct
        s = s + " -nstruct " + str(self.options.nstruct)

        #Outpath
        s = s + " -out:path:all " + self._get_make_out_path(*args, **kwargs)
        s = s + self.base_options.get_base_rosetta_flag_string()
        #Log Dir:
        if not self.options.one_file_mpi:
            dir = self._get_make_mpi_tracer_dir(*args, **kwargs)
            s = s + " -mpi_tracer_to_file "+ dir+"/tracer_logs_"

        #For these benchmarks, there are only a single root directory.
        if self.options.json_run:
            s = s + self.extra_options.get_base_rosetta_flag_string(self.base_options._get_root())
        #DB Mode
        if self.db_mode:
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

                print "Checking Features"
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

        if self.options.script_vars:
            s = s + " -parser:script_vars "+" ".join(self.options.script_vars)

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

    def _get_full_cmd(self, *args, **kwargs):
        """
        Get the full command line.
        :param args:
        :param kwargs:
        :rtype: str
        """
        cmd_string = self._get_output_string(*args, **kwargs)

        mpiexec="mpiexec"
        if self.options.mpiexec and not self.local_run():
            mpiexec = self.options.mpiexec

        if self.options.job_manager == "slurm":
            cmd = "cd "+ self._get_root() + " \n" + mpiexec + " "
        else:
            cmd = "cd "+ self._get_root() + " \n" + mpiexec + " -np " + str(self.options.np)

        if self.options.machine_file:
            cmd = cmd + " --machine_file "+self.options.machine_file+" "+ cmd_string
        else:
            cmd = cmd + " "+ cmd_string

        if self.db_mode and self.options.db_mode == "sqlite3":
            cmd = cmd + "\n"
            cmd = cmd + "cd "+self.options.outdir+"\n"
            cmd = cmd + "bash "+ get_rosetta_features_root()+"/sample_sources/merge.sh "+self.options.db_name + " "+self.options.db_name+"_*\n"
            cmd = cmd + "cd -"


        return cmd

    def local_run(self, *args, **kwargs):
        """
        Get if we are running locally
        :rtype: bool
        """
        if (self.options.job_manager == "local" or self.options.job_manager == "local_test"):
            print "Local Run!"
            return True

        else:
            return False

    def run(self, *args, **kwargs):

        log_dir = self._get_make_mpi_tracer_dir(*args, **kwargs)
        outpath = self._get_make_out_path(*args, **kwargs)
        queue_dir = self._get_make_queue_dir(*args, **kwargs)

        print "\nLogDir: "+log_dir
        print "QueueDir: "+queue_dir
        print "OutPath: "+outpath+"\n\n"


        cmd = self._get_full_cmd()

        if self.options.job_manager == "local" and self.options.print_only:
            print cmd + "\n"

        elif self.options.job_manager == "local":
            print cmd + "\n"
            os.system(cmd)
        elif self.options.job_manager == "local_test" or self.options.local_test:
            self.options.np = 2
            self.options.nstruct = 1
            self.options.split_mpi_output = False

            cmd = self._get_full_cmd()
            print_full_cmd(cmd)
            print(cmd + "\n")
            os.system(cmd)

        elif self.options.job_manager == "qsub":
            run_on_qsub(cmd, queue_dir, self._get_job_name(*args, **kwargs), self.options.print_only, self._get_job_manager_opts())

        elif self.options.job_manager == "slurm":
            run_on_slurm(cmd, queue_dir, self._get_job_name(*args, **kwargs), ntasks=self.options.np, print_only=self.options.print_only, extra_opts=self._get_job_manager_opts())


if __name__ == "__main__":
    run_rosetta = RunRosetta()
    run_rosetta.run()
