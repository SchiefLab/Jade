
A collection of mainly python modules and scripts written over the years for various purposes.

# Setup

Nothing fancy yet.  

1) Add the root path to your PYTHONPATH environment variable in your shell. 

<code>export PYTHONPATH=$PYTHONPATH:/path/to/module_c</code>

2) Add the path to module_c/bin to your PATH environment variable to use scripts and programs as executables 

<code>export PATH=$PATH:/path/to/module_c/bin</code>

# Notable Scripts and Programs

## RunRosettaMPI

Run MPI-built Rosetta locally, or an a cluster using slurm or qsub as the job manager.  Run from your root project directory or set the root dir as an option in the program.  Will cd into the root, or set the job manager script to cd into root before the MPI run.

It uses JSON files to setup the base flags (<code>--json_base</code>) and then specific flags for different rosetta runs (<code>--json_run</code>).  The [default baseline json](https://github.com/SchiefLab/module_c/blob/master/rosetta_general/jsons/common_flags.json) should be good for most runs.  See [this dir](https://github.com/SchiefLab/module_c/tree/master/rosetta_general/jsons) for a list of currently implemented jsons.  Feel free to implement your own.  I typically add that json path as an alias to easily run scripts.  The class is easily extendable for benchmarking experiments, [like I have done for antibody design](https://github.com/SchiefLab/module_c/blob/master/bin/BenchmarkRAbD.py).

Use <code>--print_only</code> to print instead of run to double check everything.  Paths can (and should) be relative.  Will setup any directories mentioned.  You can feed additional flags files or options (or overwrite any set in the json files) using :     <code>--extra_options @rel/path/to/flags rosetta_opt=setting another_opt=setting a_boolean_opt</code>

Set the job manager using the option <code>--job_manager</code>. Current options are __slurm__, __qsub__, and __local__.  Set extra options for the job manager in parenthesis, such as the slurm partition option -p, using <code>--job_manager_opts "set of -options -for run"</code>

Be sure to set <code>--np</code> and <code>--nstruct</code> (if not set in flags files or extra_options)

If you think a GUI would be useful for this, let me know!

See below for the current full help of the program:


```
usage: This program runs Rosetta MPI locally or on a cluster using slurm or qsub. in
       Relative paths are accepted.
       
       [-h] [--program PROGRAM] [--np NP] [--nodes NODES] [--ppn PPN]
       [--nstruct NSTRUCT] [-s S] [-l L] [--outdir OUTDIR]
       [--compiler {gcc,clang}] [--job_manager {slurm,qsub,local}]
       [--job_manager_opts JOB_MANAGER_OPTS] [--machine_file MACHINE_FILE]
       [--print_only] [--json_base JSON_BASE] [--json_run JSON_RUN]
       [--root ROOT] [--job_name JOB_NAME]
       [--extra_options [EXTRA_OPTIONS [EXTRA_OPTIONS ...]]] [--one_file_mpi]

optional arguments:
  -h, --help            show this help message and exit
  --program PROGRAM     Define the Rosetta program to use if not set in
                        json_run
  --np NP               
  --nodes NODES         
  --ppn PPN             Processors per node for qsub. NTasks is np for slurm
  --nstruct NSTRUCT
  -s S                  Path to a pdb file.
  -l L                  Path to a list of pdb files
  --outdir OUTDIR, -o OUTDIR
                        Outpath. Default = 'pwd/decoys'
  --compiler {gcc,clang}, -c {gcc,clang}
                        Set the compiler used. Will set clang automatically
                        for macos. Default = 'gcc'
  --job_manager {slurm,qsub,local}
                        Job Manager to launch job. Default = 'slurm '
  --job_manager_opts JOB_MANAGER_OPTS
                        Extra options for the job manager, such as queue or
                        processor requests
  --machine_file MACHINE_FILE
                        Optional machine file for passing to MPI
  --print_only          Do not actually run anything. Just print setup for
                        review.
  --json_base JSON_BASE
                        JSON file for setting up base paths/etc. for the
                        cluster.Default = 'file_dir/jsons/common_flags.json'
  --json_run JSON_RUN   JSON file for specific Rosetta run.
  --root ROOT           Set the root directory. Default is to use pwd.
                        (Benchmarking: Override any set in json_base.)
  --job_name JOB_NAME   Set the job name used for mpi_tracer_to_file dir and
                        queue. Default = 'rosetta_run'. (Benchmarking:
                        Override any set in json_base.)
  --extra_options [EXTRA_OPTIONS [EXTRA_OPTIONS ...]], -e [EXTRA_OPTIONS [EXTRA_OPTIONS ...]]
                        Extra Rosetta options. Specify like:
                        cdr_instructions=my_file other_option=setting. Note NO
                        - charactor. Booleans do not need an = sign.
  --one_file_mpi        Don't setup mpi_tracer_to_file.

```

## score_analysis

Analyze Rosetta decoys that were scored with an output json file.  Get top models, score summaries, top_n_by_10, and (soon) output pymol sessions.

Use <code>-scorefile_format json</code> during your Rosetta runs.  This is a fork of the scorefile.py script that is located in rosetta source dir.   

I copy the current help text below.

```
usage: This utility parses and extracts data from score files in JSON format
       [-h] [-s [SCORETYPES [SCORETYPES ...]]] [-n TOP_N]
       [--top_n_by_10 TOP_N_BY_10]
       [--top_n_by_10_scoretype TOP_N_BY_10_SCORETYPE] [--decoy_names] [-S]
       [--list_scoretypes]
       [scorefiles [scorefiles ...]]

positional arguments:
  scorefiles            A list of scorefiles

optional arguments:
  -h, --help            show this help message and exit
  -s [SCORETYPES [SCORETYPES ...]], --scoretypes [SCORETYPES [SCORETYPES ...]]
                        list of score terms to extract
  -n TOP_N, --top_n TOP_N
                        Only list Top N when doing top scoring decoysDefault
                        is to print all of them.
  --top_n_by_10 TOP_N_BY_10
                        Top N by 10 percent total score to print out.
  --top_n_by_10_scoretype TOP_N_BY_10_SCORETYPE
                        Scoretype to use for any top N by 10 printing. If
                        scoretype not present, won't do anything.
  --decoy_names         List decoy names
  -S, --summary         Compute stats summarizing data
  --list_scoretypes     List score term names
  
```

### Current Limitations

Works on individual scorefiles, with no -best-of-all- or combined output.


## PyRAbD_Compare

GUI for antibody design analysis.  Inputs are Antibody Features Reporter databases.  I will probably change the name soon. Each design strategy should have its own database.  Example: <code>PyRAbD_Compare.py path/to/directory/of/sqlite3/databases</code>

### Current Limitations

Note that it currently only supports sqlite3 databases and each decoy used in the comparison must have a unique name.  


# Code Organization

## _antibody_

A small collection of general antibody scripts and modules from PyIgClassify.  http://dunbrack2.fccc.edu/PyIgClassify/.  The meat of PyIgClassify should be publically released soon.


## _pymol_

Python PyMol modules and pymol scripts


## _rosetta_general_

Rosetta (www.rosettacommons.org) modules and flags files for analyzing results, benchmarking, etc.  PyRosetta (www.pyrosetta.org) modules and scripts from various projects


## _sequence_

Modules for dealing with protein sequence


## _structure_

Modules for reading PDBs and storing structure information.  Yes, my own general PDB reader.  Because everyone has one, right?


## _database_

Text files, jsons, etc. for import into scripts and modules.


## _tcl_

TCL modules for molecular dynamic simulations.


