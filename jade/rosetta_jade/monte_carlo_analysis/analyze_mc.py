from collections import defaultdict
import sys, os, re
import pandas
from basic import path

# Utility file for calculating mc data. Currently, we have no good c++ methods to do this

def read_mc_acceptance_from_pdbs(pdb_paths, list_of_pdbids = []):
    """
    Reads acceptance from a list of PDBs, returns a pandas dataframe for further processing.

    pdbid list is optional to create a new column in the dataframe.  This is useful for benchmarking.
    It will match the pdbid with the filename, adding it as a new column.

    :param pdb_paths: list
    :param list_of_pdbids: list
    :rtype: pandas.DataFrame

    """
    for pdb_path in pdb_paths:
        if not os.path.exists(pdb_path):
            sys.exit("Path does not exist: "+pdb_path)

        INFILE = path.open_file(pdb_path, 'r')

        for line in INFILE:
            line = line.strip()
            if not line or line.startswith('#'):continue


            pdb_id = get_pdb_id_of_path(file_path)
            lineSP = line.split()

            if not lineSP:                continue
            if not lineSP[0] == "ACCEPT": continue
            if not lineSP[1] == "LOG":    continue

            if lineSP[3] == "FINAL":
                if lineSP[4] == "END":
                  cy_data.set_protocol_final_e(lineSP[5])
                else:
                  cy_data.set_final_energy(lineSP[4], lineSP[5])
            else:
             if lineSP[3] == "0":
                  cy_data.set_protocol_start_e(lineSP[4])
             elif lineSP[3] == "-1":
                   cy_data.set_protocol_native_e(lineSP[4])
              else:
                  cy_data.set_mc_data(lineSP[3], lineSP[4], lineSP[5])

        INFILE.close()

