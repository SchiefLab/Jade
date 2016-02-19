#!/usr/bin/env python

#Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import os
import gzip
import re

from Bio.PDB import PDBParser

from Bio.PDB import MMCIFParser

from basic.RestypeDefinitions import RestypeDefinitions
from basic.structure.Structure import PDBResInfo as PDBInfo
from basic.structure.Structure import Residue as res_struct
from utility import vector1


class BioPose(object):
    """
    This is my biopython meta class.  Because biopython's interface kinda sucks.
    This is a little cleaner.

    The other way is to sublclass each Biopython class structure, which I'm not ready to do.

    Right now, you need a path as I don't know how we would use this from sequence, etc as you do in Rosetta.
    """
    def __init__(self, path):

        self.struct, self.header = self.load_from_file(path) #Bio struct, Header dictionary

        self.all_residues = self._setup_all_residues(model_num=0) #vector1 of Bio Residues
        self.pdb_info = self._setup_pdb_info(model_num=0) #PDBInfo to map the vector1

        self.res_definitions = RestypeDefinitions()

    ############ IO ###################
    def load_from_file(self, path):
        """
        Load a file from PDB or mmCIF.  .gz is supported.

        :param path: Path to PDB or mmCIF file
        :rtype: tuple(bio.PDB.Structure.Structure, dict)
        """
        structure = None
        if re.search(".pdb", path):
            parser = PDBParser()
        else:
            parser = MMCIFParser()

        path = path.strip()
        model_id = os.path.basename(path)
        if os.path.basename(path).split('.')[-1] == 'gz':
            GZ = gzip.open(path, 'rb')
            GZ.close()
        else :
            structure = parser.get_structure(model_id, path)
        header = parser.get_header()

        return structure, header

    def reload_from_file(self, path):
        self.struct, self.header = self.load_from_file(path)

    ############ Getting Components #############
    def structure(self):
        """
        Get the Bio Structure stored in this class.
        :rtype: bio.PDB.Structure.Structure
        """
        return self.struct

    def model(self, model_num = 0):
        """
        Get a Bio Model of the stored structure
        :param id: int
        :rtype: bio.PDB.Model.Model
        """
        return self.struct[model_num]

    def chain(self, chain_id, model_num = 0):
        """
        Get a Bio Chain of the stored structure
        :param chain_id: str
        :param model_num: int
        :rtype: bio.PDB.Chain.Chain
        """
        return self.struct[model_num][chain_id]

    def residue(self, resnum, chain_id, icode=' ', alt=' ', model_num = 0):
        """
        Get a Bio Residue of the stored structure

        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: bio.PDB.Residue.Residue
        """
        return self.struct[model_num][chain_id][(alt, resnum, icode)]


    def atom(self, atom_name, resnum, chain_id, icode=' ', alt=' ', model_num=0):
        """
        Get a Bio Atom of the stored structure

        :param atom_name: str
        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: bio.PDB.Atom.Atom
        """
        return self.struct[model_num][chain_id][(alt, resnum, icode)][atom_name]


    ### Lists of Structures ###

    def chains(self, model_num = 0):
        """
        Get a list of Bio Chains
        :param model_num: int
        :rtype: list[bio.PDB.Chain.Chain
        """
        return [c for c in self.struct]

    def residues(self, chain_id, model_num = 0, include_alt = False):
        """
        Get residues, including or not including residues with alternate location codes - which can be a PITA

        :param chain_id: str
        :param model_num: int
        :param include_alt: bool
        :rtype: list[bio.PDB.Residue.Residue]
        """
        resi = []
        for res in self.chain(chain_id, model_num):
            if res.id[0] ==' ':
                resi.append(res)
            elif include_alt:
                resi.append(res)
            else:
                continue
        return resi

    def atoms(self, resnum, chain_id, icode=' ', alt=' ', model_num = 0):
        """
        Get a list of Bio Atoms
        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: list[bio.PDB.Atom.Atom]
        """
        return [atm for atm in self.residue(resnum, chain_id, icode, alt, model_num)]


    ############ Helper Funtions ################
    def get_sequence(self, chain_id, model_num = 0):
        """
        Get a sequence of a chain - Not including alternate res locations

        :param chain_id: str
        :param model_num: int
        :rtype: str
        """
        if self.get_chain_length(chain_id, model_num) == 0:
            return ""

        seq = ""
        for res in self.residues(chain_id, model_num):
            aa = self.res_definitions.get_one_letter_from_three(res.resname)
            if not aa:
                print "Setting NCAA as X: "+res.resname
                print "This could pose a problem!"
                seq = seq+'X'
                continue

            seq = seq+aa
        return seq

    def get_chain_length(self, chain_id, model_num = 0):
        """
        Get the number of AA in a chain - Not including alternate res locations
        :param chain_id: str
        :rtype: int
        """
        return len(self.residues(chain_id, model_num))

    def get_chain_ids(self, model_num):
        """
        Get all chain IDS for a model.
        :param model_num: int
        :rtype: list[str]
        """
        ids = []
        for chain in self.model(model_num):
            ids.append(chain.id)
        return ids


    ###### Private Members #######

    def _setup_all_residues(self, model_num=0):
        """
        Setup all residues in the model, for ease of use like a pose.
        :param model_num: int
        :rtype: list[bio.PDB.Residue.Residue]
        """
        all_residues = vector1()

        for chain_id in self.get_chain_ids(model_num):
            residues = self.residues(chain_id, model_num)
            all_residues.extend(residues)

        return all_residues

    def _setup_pdb_info(self, model_num=0):
        """
        Setup the PDBInfo mapping from residues to pdb iformation.
        :param model_num: int
        :rtype: PDBInfo
        """
        pdb_info = PDBInfo()
        i = 1
        for res in self.all_residues:
            s = res_struct(self.res_definitions.get_one_letter_from_three(res.resname), res.id[1], res.id[2])
            pdb_info.add_residue(s)

        return pdb_info


if __name__ == "__main__":
    v = vector1([1, 2, 3])
    for i in v:
        print repr(i)

    print len(v)

    print repr(v[3])






