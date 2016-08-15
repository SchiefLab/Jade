

A repository for modules and applications to aid in the design and analysis of Biological molecules, especially when working with Rosetta or PyRosetta.

# Authors

Jared Adolf-Bryfogle (jadolfbr@gmail.com), lab of Dr. William Schief

# Contributions

If you have great changes you would like to add to the repository, please open a Pull Request.  Once people actually start using this, I will switch to this model myself via branches.

Currently, we do not have any automatic testing in place.  Unit testing is planned.  Please be careful with what you code.  If you are intrepid, you can go ahead and add unit testing to /testing.  This will be done when I have time.

# License 

Like Numpy, and SciPy (amongst many other python modules), we use a [BSD license](https://opensource.org/licenses/BSD-2-Clause).


# Setup

Nothing fancy yet.  A true python install via PIP is planned.  For now, you will need to set it up manually. 

1) __PYTHONPATH:__ Add the root path to your _PYTHONPATH_ environment variable in your shell. 

 - <code>export PYTHONPATH=$PYTHONPATH:/path/to/Jade/src</code>

2) __PATH:__ Add the path to Jade/apps to your _PATH_ environment variable to use scripts and programs as executables 

 - <code>export PATH=$PATH:/path/to/Jade/apps</code>

3) __Dependancies:__  
 - Install [Pip](https://pypi.python.org/pypi/pip) if you don't already have it.  
 - Run: ```sudo ./setup_dependancies.sh``` in the main Jade directory.
 
4) __OptionalL:__ 
 - _RABD (RosettaAntibodyDesign) Applications_: Add the path to Jade/apps/RAbD to your PATH environment variable.
 
  -  ```export PATH=$PATH:/path/to/Jade/apps/RAbD```
 
 
 
# Code Organization

<code>Jade/apps</code>
 - Applications, and scripts
  
<code>Jade/database</code>
 - Collection of files used by Jade applications and modules.
 
<code>Jade/src</code>
 - Jade Source Code
 
<code>Jade/testing</code>
 - Testing code and inputs.  Not yet developed fully.





# Jade SRC Code

## _basic_
Useful general classes and collections of functions (Threading, BioPose, PandasDataFrame, path, etc)

## _utility_
Functions and simple classes go in <code>__init__.py</code> 
vector1 is a list indexed at 1

 - Use: <code>from utility import vector1</code>

## _antibody_
A small collection of general antibody scripts and modules from PyIgClassify.  http://dunbrack2.fccc.edu/PyIgClassify/.  The meat of PyIgClassify should be publically released soon.

## _plotting_
Collection of plotting classes and functions for matplotlib, seaborn

## _pymol_jade_
Python PyMol modules and pymol scripts


## _rosetta_jade_
Rosetta (www.rosettacommons.org) modules and flags files for analyzing results, benchmarking, etc.  PyRosetta (www.pyrosetta.org) modules and scripts from various projects


## _sequence_
Modules for dealing with protein sequence


## _structure_
Modules for reading PDBs and storing structure information.  Yes, my own general PDB reader.  Because everyone has one, right?

## _tcl_
TCL modules for molecular dynamic simulations.


# Notable Scripts and Programs

## score_analysis

Analyze Rosetta decoys that were scored with an output json file.  Get top models, score summaries, top_n_by_10, and output pymol sessions, plot scores.

Use <code>-scorefile_format json</code> during your Rosetta runs (but now currently works on json or classic score files.  This is a fork of the scorefile.py script that is located in rosetta source dir.   

I copy the current help text below.

```
usage: This utility parses and extracts data from score files in JSON format
       [-h] [-s [SCORETYPES [SCORETYPES ...]]] [-n TOP_N]
       [--top_n_by_10 TOP_N_BY_10]
       [--top_n_by_10_scoretype TOP_N_BY_10_SCORETYPE]
       [--decoy_names [DECOY_NAMES [DECOY_NAMES ...]]] [--list_scoretypes]
       [--pdb_dir PDB_DIR] [--summary] [--csv] [--make_pdblist]
       [--pymol_session] [--plot [PLOT [PLOT ...]]] [--copy_top_models]
       [--prefix PREFIX] [--outdir OUTDIR]
       [--plot_type {line,scatter,bar,hist,box,kde,area,pie,hexbin}]
       [--plot_filter PLOT_FILTER] [--native NATIVE] [--ab_structure]
       [--super SUPER]
       [scorefiles [scorefiles ...]]

positional arguments:
  scorefiles            A list of scorefiles

optional arguments:
  -h, --help            show this help message and exit
  -s [SCORETYPES [SCORETYPES ...]], --scoretypes [SCORETYPES [SCORETYPES ...]]
                        List of score terms to extract
  -n TOP_N, --top_n TOP_N
                        Only list Top N when doing top scoring decoys or
                        making pymol sessionsDefault is to print all of them.
  --top_n_by_10 TOP_N_BY_10
                        Top N by 10 percent total score to print out.
  --top_n_by_10_scoretype TOP_N_BY_10_SCORETYPE
                        Scoretype to use for any top N by 10 printing. If
                        scoretype not present, won't do anything.
  --decoy_names [DECOY_NAMES [DECOY_NAMES ...]]
                        Decoy names to use
  --list_scoretypes     List score term names
  --pdb_dir PDB_DIR, -d PDB_DIR
                        Directory for PDBs if different than the directory of
                        the scorefile

OUTPUT:
  General output options.

  --summary, -S         Compute stats summarizing data
  --csv, -c             Output selected columns, top, and decoys as CSV.
  --make_pdblist        Output PDBlist file(s)
  --pymol_session       Make pymol session(s) of the scoretypes specified
  --plot [PLOT [PLOT ...]]
                        Plot one score type vs another. Save the plot. 2 or 3
                        Arguments. [X, Y, 'Title''] OR [X, 'Title']. If title
                        has spaces, use quotes. Nothing special, just used for
                        quick info.
  --copy_top_models     Copy the top -n to the output directory for each
                        scorefile passed.
  --prefix PREFIX, -p PREFIX
                        Prefix to use for any file output. Do not include any
                        _
  --outdir OUTDIR, -o OUTDIR
                        Output dir. Default is current directory.

PLOTTING:
  Options for plot output

  --plot_type {line,scatter,bar,hist,box,kde,area,pie,hexbin}
                        The type of plot we are outputting.
  --plot_filter PLOT_FILTER
                        Filter X to top Percent of this - useful to remove
                        outliers.

PYMOL:
  Options for pymol session output

  --native NATIVE       Native structure to use for pymol sessions.
  --ab_structure        Specify if the module is a renumbered antibody
                        structure. Will run pymol script for ab-specific
                        selection
  --super SUPER         Super this selection instead of align all to.
  
```
### Current Limitations

Works on individual scorefiles, with no -best-of-all- or combined output.

## RunRosettaMPI

Run MPI-built Rosetta locally, or an a cluster using slurm or qsub as the job manager.  Run from your root project directory or set the root dir as an option in the program.  Will cd into the root, or set the job manager script to cd into root before the MPI run.

It uses JSON files to setup the base flags (<code>--json_base</code>) and then specific flags for different rosetta runs (<code>--json_run</code>).  The [default baseline json](https://github.com/SchiefLab/module_c/blob/master/rosetta_jadeeral/jsons/common_flags.json) should be good for most runs.  See [this dir](https://github.com/SchiefLab/module_c/tree/master/rosetta_jadeeral/jsons) for a list of currently implemented jsons.  Feel free to implement your own.  I typically add that json path as an alias to easily run scripts.  The class is easily extendable for benchmarking experiments, [like I have done for antibody design](https://github.com/SchiefLab/module_c/blob/master/bin/BenchmarkRAbD.py).

Use <code>--print_only</code> to print instead of run to double check everything.  Paths can (and should) be relative.  Will setup any directories mentioned.  You can feed additional flags files or options (or overwrite any set in the json files) using :     <code>--extra_options @rel/path/to/flags rosetta_opt=setting another_opt=setting a_boolean_opt</code>

Set the job manager using the option <code>--job_manager</code>. Current options are __slurm__, __qsub__, and __local__.  Set extra options for the job manager in parenthesis, such as the slurm partition option -p, using <code>--job_manager_opts "set of -options -for run"</code>

Be sure to set <code>--np</code> and <code>--nstruct</code> (if not set in flags files or extra_options)

Relational Database support has been added.  Use RunRosettaDBMode app instead.  If using sqlite3, it will automatically combine the databases at the end of the run.  Very useful for running features reporters.  

If you think a GUI would be useful for this, let me know!
See below for the current full help of the program:


```
usage: This program runs Rosetta MPI locally or on a cluster using slurm or qsub.  Relative paths are accepted.
       [-h] [-s S] [-l L] [--np NP] [--nstruct NSTRUCT] [--job_name JOB_NAME]
       [--outdir OUTDIR] [--json_run JSON_RUN] [--extra_options EXTRA_OPTIONS]
       [--script_vars [SCRIPT_VARS [SCRIPT_VARS ...]]] [--program PROGRAM]
       [--print_only] [--local_test] [--one_file_mpi]
       [--job_manager {slurm,qsub,local,local_test}]
       [--job_manager_opts [JOB_MANAGER_OPTS [JOB_MANAGER_OPTS ...]]]
       [--json_base JSON_BASE] [--compiler {gcc,clang}] [--mpiexec MPIEXEC]
       [--machine_file MACHINE_FILE]

optional arguments:
  -h, --help            show this help message and exit

Common Options:
  -s S                  Path to a pdb file
  -l L                  Path to a list of pdb files
  --np NP               Number of processors to use for MPI. Default = 101
  --nstruct NSTRUCT     The number of structures/parallel runs. Can also set
                        this in any JSON file.
  --job_name JOB_NAME   Set the job name used for mpi_tracer_to_file dir and
                        queue. Default = 'rosetta_run'. (Benchmarking:
                        Override any set in json_base.)
  --outdir OUTDIR, -o OUTDIR
                        Outpath. Default = 'pwd/decoys'
  --json_run JSON_RUN   JSON file for specific Rosetta run. Not required. Pre-
                        Configured JSONS include: ['antibody_designer.json',
                        'antibody_designer_dock.json',
                        'antibody_designer_even_clus.json',
                        'antibody_designer_even_clus_dock.json',
                        'antibody_designer_even_len_clus.json',
                        'antibody_designer_even_len_clus_dock.json',
                        'antibody_features.json', 'blank.json',
                        'cluster_features.json', 'common_flags.json',
                        'dualspace_relax.json', 'glycan_clash_check.json',
                        'glycosylate_relax.json', 'interface_analyzer.json',
                        'NGK.json', 'NGK_smooth.json', 'NGK_smooth_shap.json',
                        'pareto_optimal_relax.json', 'relax.json',
                        'relaxed_design.json', 'relaxed_design_ds.json',
                        'remodel.json', 'rosetta_scripts.json']
  --extra_options EXTRA_OPTIONS
                        Extra Rosetta options. Specify in quotes!
  --script_vars [SCRIPT_VARS [SCRIPT_VARS ...]]
                        Any script vars for XML scripts.Specify as you would
                        in Rosetta. like: glycosylation=137A,136A
  --program PROGRAM     Define the Rosetta program to use if not set in
                        json_run

Testing and Debugging:
  --print_only          Do not actually run anything. Just print setup for
                        review.
  --local_test          Is this a local test? Will change nstruct to 1 and run
                        on 2 processors
  --one_file_mpi        Output all MPI std::out to a single file instead of
                        splitting it.

Special Options for controlling execution:
  --job_manager {slurm,qsub,local,local_test}
                        Job Manager to launch job. (Or none if local or
                        local_test)Default = 'slurm '
  --job_manager_opts [JOB_MANAGER_OPTS [JOB_MANAGER_OPTS ...]]
                        Extra options for the job manager, such as queue or
                        processor requestsRemove double dashes. Exclusive is
                        on by default. Specify like: -p imperial exclusive.
  --json_base JSON_BASE
                        JSON file for setting up base paths/etc. for the
                        cluster.Default =
                        'database/rosetta/jsons/common_flags.json'
  --compiler {gcc,clang}
                        Set the compiler used. Will set clang automatically
                        for macos. Default = 'gcc'
  --mpiexec MPIEXEC     Specify a particular path (or type of) MPI exec.
                        Default is srun (due to vax). If local or local test,
                        will use mpiexex
  --machine_file MACHINE_FILE
                        Optional machine file for passing to MPI

```


## get_seq

My own app to get sequences from structures for regions, chains, etc, and output fasta files or just print to the screen (either as a single PDB or multiple).  Also used to create Gene order forms for genscript.  THis is a bit specific for the SchiefLab, but it can be generalized if needed. 

```
usage: get_seq.py [-h] [--pdb PDB] [--pdblist PDBLIST]
                  [--pdblist_input_dir PDBLIST_INPUT_DIR] [--chain CHAIN]
                  [--cdr CDR]
                  [--format {basic,fasta,general_order,IgG_order,IgG_order_lambda,IgG_order_kappa,IgG_order_heavy}]
                  [--outpath OUTPATH] [--prefix PREFIX] [--region REGION]
                  [--strip_c_term STRIP_C_TERM] [--pad_c_term PAD_C_TERM]
                  [--output_original_seq]

Uses Biopython to print sequence information. Example: get_seq.py --pdb
2j88_A.pdb --format fasta --outpath test.txt

optional arguments:
  -h, --help            show this help message and exit
  --pdb PDB, -s PDB     Input PDB path
  --pdblist PDBLIST, -l PDBLIST
                        Input PDB List
  --pdblist_input_dir PDBLIST_INPUT_DIR, -i PDBLIST_INPUT_DIR
                        Input directory if needed for PDB list
  --chain CHAIN, -c CHAIN
                        A specific chain to output
  --cdr CDR             Pass a specific CDR to output alignments of.
  --format {basic,fasta,general_order,IgG_order,IgG_order_lambda,IgG_order_kappa,IgG_order_heavy}
                        The output format requried.
  --outpath OUTPATH, -o OUTPATH
                        Output path. If none is specified it will write to
                        screen.
  --prefix PREFIX, -t PREFIX
                        Tag to add before chain
  --region REGION       specify a particular region, start:end:chain
  --strip_c_term STRIP_C_TERM
                        Strip this sequence off the C-term of resulting
                        sequences. (Useful for antibodies
  --pad_c_term PAD_C_TERM
                        Pad this sequence with some C-term (Useful for
                        antibodies
  --output_original_seq
                        Output the original sequence and the striped seqeunce
                        if stripped. Default FALSE.
 ```

