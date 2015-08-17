#Jared Adolf-Bryfogle
#Reorders PDBFiles in a dirctory according to LH_A in order for Rosetta Antibody Design benchmarking. Removes HetAtm!!!

import os
import sys

from structure.PythonPDB2 import *

def reorder_and_save_chains(in_path, out_path, remove_het = False):
    blank_pdb = PythonPDB2()
    full_pdb = PythonPDB2(in_path)

    blank_pdb.copy_chain_into_pdb_map(full_pdb, "A")
    blank_pdb.copy_chain_into_pdb_map(full_pdb, "L")
    blank_pdb.copy_all_but_chains_into_pdb_map(full_pdb, ["A", "L"])

    if remove_het:
        blank_pdb.remove_hetatm_atoms()
        blank_pdb.remove_waters()

    blank_pdb.save_PDB(out_path)

if __name__ == "__main__":

    #I don't have time to make this fancy right now.

    #Arguments:
    ### Input Directory
    ### Input PDBList
    ### Output Directory

    #Assumes PDBList has no paths and requires an input directory as if we run Rosetta.

    in_dir = sys.argv[1]
    in_pdblist = sys.argv[2]
    out_dir = sys.argv[3]

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    PDBLIST = open(in_pdblist, 'r')
    for line in PDBLIST:
        line = line.strip()
        line = in_dir+"/"+line

        print "Reordering "+line
        outpath = out_dir+"/"+os.path.basename(line)
        reorder_and_save_chains(line, outpath)

    print "Done!!"
    PDBLIST.close()