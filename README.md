
A collection of mainly python modules and scripts written over the years for various purposes.

# Use

Nothing fancy yet.  Just add the root path to your PythonPath in your shell. 

# Useful Modules and Scripts

## _antibody_

A small collection of general antibody scripts and modules from PyIgClassify.  http://dunbrack2.fccc.edu/PyIgClassify/.  The meat of PyIgClassify should be publically released soon.

* __ab_db.py__

   *   Functions for interacting with the PyIgClassify antibody database.

* __split_antibody_components.py__

   *   Split antibody into Fv, Fc, chains, etc. through the use of AHo numbering


## _pymol_

Python PyMol modules and pymol scripts

* __color_cdrs.pml__

   *   Color CDRs based on AHo numbering.  I usually copy this script to a new root directory at /pymol and then use in pymol by @/pymol/color_cdrs.pml
   
* __PyMolScriptWriter.py__
   
   *  Module for writing out loads of pymol scripts and gerating pymol sessions through pymol code using this class.  
   
   
## _rosetta_

Rosetta (www.rosettacommons.org) modules and flags files for analyzing results, benchmarking, etc.  PyRosetta (www.pyrosetta.org) modules and scripts from various projects

* __alignment.py__
   *  A collection of PyRosetta alignment functions
   
* __DesignBreakdown.py__
   *  Output protein design sequence statistics in text and database form.  Plot in R.
   
* __get_mutation_energy.py__
   *  Scan aa mutations for a set of positions, relaxing at each one or not.  Or do an alanine scan using talaris2014 energy function.  PyRosetta
   
* __RunRosetta.py__ and __SetupRosettaOptionsGeneral.py__
   *  Classes for setting up and running rosetta on a cluster using MPI and dealing with directories, pathing, and queue management.   Used for benchmarking via subclass. 


## _sequence_

Modules for dealing with protein sequence

* __ClustalRunner.py__
   * Class wrapper for running the Clustal program
   
* __PDBConsensusInfo.py__ and __SequenceStats.py__
   * Basic module for getting and storing frequency and probability info and writing to sequence logos

## _structure_

Modules for reading PDBs and storing structure information.  Yes, my own general PDB reader.  Because everyone has one, right?

* __PythonPDB2.py__
   *  My PDB reading module.  Light and simple.
  
* __SQLPDB.py__
   *  Read and write PDB information to/from sqlite3 databases

* __Structure.py__
   *  Representations of protein and molecular structure.  Mainly just for holding and accessing information

## _root_

* __restype_definitions.py__
   *  Amino acid sequence information, different names, name1, name3, full names, etc. A way to do this outside Rosetta.
  
* __calibur.py__
   *  A class wrapper for running the calibur clustering application.
   
* __fasta.py__
   * A collection of functions for dealing with sequence information and fast files.  Requires PyRosetta and Biopython



## _database_

Text files, jsons, etc. for import into other programs




## _tcl_

TCL modules for molecular dynamic simulations.


