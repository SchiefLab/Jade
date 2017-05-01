#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/Structure.py
## @brief  Class for structural representations of specific protein types.  Antibody and CDR work, feel free to add.  CDR Stuff Will be in C++ Rosetta soon.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys
from collections import defaultdict
from jade.basic.RestypeDefinitions import RestypeDefinitions

### Global Values #####
cdr_names = ["L1", "L2", "L3", "H1", "H2", "H3"]
heavy_cdr_names = ["H1", "H2", "H3"]



#########################################
class PDBInfo(object):
    """
    Analogous to Rosetta PDBInfo Class
    I should start at 1
    """

    def __init__(self):
        self.pose_to_record_map = defaultdict()
        self.pdb_to_pose_map = defaultdict()
        self.extra_data = defaultdict()

    def clear(self):
        self.pose_to_record_map = defaultdict()
        self.pdb_to_pose_map.clear()
        self.extra_data = defaultdict()

    def set_residue_record(self, i, residue_record):
        if not isinstance(residue_record, ResidueRecord):
            sys.exit("Residue class must be passed for residue_record!")
        self.pose_to_record_map[i] = residue_record
        self.pdb_to_pose_map[ (residue_record.tuple())] = i


    def add_residue_record(self, residue_record):
        if not isinstance(residue_record, ResidueRecord):
            sys.exit("Residue class must be passed for residue_record!")
        self.pose_to_record_map[len(self.pose_to_record_map)+1] = residue_record
        self.pdb_to_pose_map[ (residue_record.tuple() )] = len(self.pose_to_record_map)

    def set_pdb_num(self, i, pdb_num, icode =" "):

        self.pdb_to_pose_map.pop(self.pose_to_record_map[i].get_tuple())

        self.pose_to_record_map[i].set_pdb_num(pdb_num)
        self.pose_to_record_map[i].set_icode(icode)

        self.pdb_to_pose_map[ self.pose_to_record_map[i].get_tuple() ] = i

    def set_icode(self, i, icode):
        self.pdb_to_pose_map.pop(self.pose_to_record_map[i].get_tuple())
        self.pose_to_record_map[i].set_icode(icode)
        self.pdb_to_pose_map[ self.pose_to_record_map[i].get_tuple() ] = i


    def get_all_residue_records(self):
        """
        return all residue_records held as an array in order.
        """
        residue_records = []
        for i in range(1, self.total_residue_record()+1):
            residue_records.append(self.get_residue_record(i))

        return residue_records

    def set_extra_data(self, key, value):
        self.extra_data[key] = value

    def get_extra_data(self, key):
        return self.extra_data[key]


    #########################
    def pose2pdb(self, i):
        rec = self.get_residue_record(i)
        return [str(rec.pdb_num), rec.chain, rec.icode]

    def pdb2pose(self, resnum, chain_id, icode=' '):
        return self.pdb_to_pose_map[(resnum, chain_id, icode)]


    #########################
    def delete_residue_record(self):
        """
        Delete the residue_record and renumber starting from 1
        """

        pass

    def get_residue_record(self, i):
        """
        Get the residue record class for this particular index.
        :param i:
        :rtype: ResidueRecord
        """
        return self.pose_to_record_map[i]

    def res(self, i):
        return self.pose_to_record_map[i]

    def get_resnum(self, pdb_num, chain, icode = ' '):
        """
        Get the matching 'resnum' (i) or None if not found.
        """
        return self.pdb2pose(pdb_num, chain, icode)

    def get_residue_record_of_pdb_num(self, pdb_num, chain_id, icode =' '):
        for residue_record in self.pose_to_record_map:
            if residue_record.pdb_num == pdb_num and residue_record.icode == icode:
                return residue_record


    def get_sequence(self):
        seq = ""
        for i in range(1, self.total_residue_record() +1):
            seq = seq + self.pose_to_record_map[ i ].get_aa()
        return seq

    def get_sequence_bt_resnums(self, start, stop):
        seq = ""
        for i in range(start, stop +1):
            seq = seq + self.pose_to_record_map[ i ].get_aa()

        print seq
        return seq


    def get_sequence_bt_residue_records(self, res1, res2, chain):
        #print repr(res1)
        #print repr(res2)

        seq = ""
        for i in range(1, self.total_residue_records()+1):
            #print repr(i)
            res = self.pose_to_record_map[ i ]
            if res.get_pdb_num() >= res1.get_pdb_num() and res.get_pdb_num() <= res2.get_pdb_num() and res.get_chain() == chain:
                seq = seq+res.get_aa()

        print seq
        return seq

    def total_residue_records(self):
        return len(self.pose_to_record_map)

    def total_residue_record(self):
        return len(self.pose_to_record_map)

    def total_residue(self):
        return self.total_residue_record()

class ResidueRecord(object):
    """
    Basic class to PDBInfo
    """
    def __init__(self, one_letter_aa, pdb_num, chain, icode = " "):
        self.aa = one_letter_aa
        self.pdb_num = pdb_num
        self.chain = chain
        self.icode = icode
        self.extra_info = defaultdict()

    def __str__(self):
        s = repr(self.pdb_num)+" "+self.chain+" "+self.icode+" "+ self.aa
        return s

    def __repr__(self):
        return self.__str__()

    def tuple(self):
        return int(self.pdb_num), self.chain, self.icode

    def set_pdb_num(self, pdb_num):
        self.pdb_num = pdb_num

    def set_chain(self, chain):
        self.chain = chain

    def set_icode(self, icode):
        self.icode = icode

    def set_aa(self, one_letter_aa):
        self.aa = one_letter_aa

    def set_extra_info(self, key, value):
        self.extra_info[key] = value

    ######################
    def get_aa(self):
        return self.aa

    def get_pdb_num(self):
        return self.pdb_num

    def get_chain(self):
        return self.chain

    def get_icode(self):
        return self.icode

    def has_extra_info(self, key):
        if self.extra_info.has_key(key):
            return True
        else:
            return False

    def get_extra_info(self, key):
        return self.extra_info[key]

    def get_extra_info_keys(self):
        return sorted(self.extra_info.keys())

    def get_extra_info_dict(self):
        return self.extra_info

    def init_extra_infos(self, array_of_keys, value):
        for key in array_of_keys:
            self.extra_info[key] = value

class AntibodyResidueRecord(ResidueRecord):
    """
    Extension of Residue used to hold and access extra data used for renumbering/printing renumbering info
    I could backport python Enums, which would be incredibly useful here, but I don't want to require the additional step.
     - used in Python3.4
    """
    def __init__(self, aa, pdb_res_num, chain, icode=" "):
        ResidueRecord.__init__(self, aa, pdb_res_num, chain, icode)
        extra_infos = ["cdr_type",
                       "old_resnum",
                       "old_chain",
                       "old_icode",
                       "region",
                       "chain_type",
                       "cluster",
                       "cluster_dis",
                       "gene",
                       "meta"]

        self.init_extra_infos(extra_infos, "")

    def is_cdr(self):
        if self.get_region() == "CDR":
            return True
        else:
            return False

    def is_framework(self):
        if self.get_region() == "FRAME":
            return True
        else:
            return False


    ### Old Numbers ###
    ###
    def set_old_resnum(self, old_resnum, icode):
        """
        Sets the old resnum and icode in the numbering map dictionary.
        """
        self.set_extra_info('old_resnum', old_resnum)
        self.set_extra_info('old_icode', icode)
    def get_old_resnum(self):
        return self.get_extra_info('old_resnum')

    def set_old_icode(self, old_icode):
        self.set_extra_info("old_icode", old_icode)
    def get_old_icode(self):
        return self.get_extra_info('old_icode')

    def set_old_chain(self, old_chain):
        self.set_extra_info('old_chain', old_chain)

    def get_old_chain(self):
        return self.get_extra_info('old_chain')

    ### Regions ###
    ###
    def set_region(self, region):
        self.set_extra_info("region", region)
    def get_region(self):
        """
        Get the region type of the position - CDR/FRAMEWORK
        """
        return self.get_extra_info('region')

    def set_chain_type(self, chain_type):
        self.set_extra_info("chain_type", chain_type)
    def get_chain_type(self):
        """
        Gets the chaintype for the position - L or H
        """
        return self.get_extra_info('chain_type')

    def set_gene(self, gene):
        self.set_extra_info("gene", gene)
    def get_gene(self):
        return self.get_extra_info('gene')

    ### CDRs ###
    ###
    def set_cdr_type(self, cdr_type):
        self.set_extra_info("cdr_type", cdr_type)

    def get_cdr_type(self):
        return self.get_extra_info('cdr_type')

    def set_cluster(self, cluster):
        self.set_extra_info("cluster", cluster)
    def get_cluster(self):
        return self.get_extra_info('cluster')

    def set_distance(self, distance):
        self.set_extra_info("distance", distance)
    def get_distance(self):
        return self.get_extra_info('cluster_dis')

    def set_meta(self, meta):
        self.set_extra_info("meta", meta)
    def get_meta(self):
        return self.get_extra_info("meta")

class ResidueRegion:
    def __init__(self, res1, res2, name = None):
        assert isinstance(res1, ResidueRecord)
        assert isinstance(res2, ResidueRecord)

        self.res1 = res1
        self.res2 = res2
        self.name = name

    def __str__(self):
        return repr(self.res1)+" : "+repr(self.res2)

    def __repr__(self):
        return self.__str__()

    def get_name(self):
        return self.name

    def get_res1(self):
        return self.res1

    def get_res2(self):
        return self.res2

    ###Lots of other cool stuff to do with this:

class AntibodyStructure:
    """
    Simple class for accessing Modified_AHO antibody numbering information outside of Rosetta.
    """
    def __init__(self):
        #3 Different ways to access the CDR Loops.  Should be one.  Not convoluted.

        self.cdrs = defaultdict()
        self.frameworks = defaultdict()
        self.cdr_names = cdr_names #So we have an ordered list the way we like it.

        for cdr_name in self.cdr_names:
            self.cdrs[cdr_name] = CDR(cdr_name)

        for chain in ["L", "H"]:
            self.frameworks[chain] = FrameworkRegions(chain)

    def get_cdr(self, cdr_name):
        """
        Get a CDR object of the given name
        :param cdr_name: str
        :rtype: CDR
        """
        return self.cdrs[cdr_name]

    def get_cdrs(self):
        """
        Get a dictionary of CDR objects, indexed by name
        :rtype: defaultdict[str] = CDR

        """
        return self.cdrs


    def get_cdr_names(self):
        """
        Get a list of cdr names
        :rtype: [str]
        """
        return self.cdr_names

    def get_cdr_seq(self, bio_pose, cdr_name):
        """
        Get the CDR sequence from a bio pose
        :param bio_pose: basic.structure.BioPose
        :param cdr_name: str
        :rtype: str
        """
        seq = ""

        start = self.get_cdr(cdr_name).get_pdb_start()
        end = self.get_cdr(cdr_name).get_pdb_end()
        chain = self.get_cdr(cdr_name).get_pdb_chain()

        return bio_pose.pdbinfo().get_sequence_bt_resnums(bio_pose.resnum(start, chain), bio_pose.resnum(end, chain))


    def get_framework_info(self, chain):
        """
        Get the framework info class
        :param chain: str
        :rtype: FrameworkRegions
        """
        return self.frameworks[chain]

class CDR:
    def __init__(self, name):
        self.regions = {
            "L1":['L', 24, 42],
            "L2":['L', 57, 72],
            "L3":['L', 107, 138],
            "H1":['H', 24, 42],
            "H2":['H', 57, 69],
            "H3":['H', 107, 138]}

        self.name = name
        self.region = self.regions[self.name]
        self.chain = self.region[0]
        self.Nter = self.region[1]
        self.Cter = self.region[2]
        self.residue_records = dict()

    def __str__(self):
        return str(self.regions[self.name])

    def get_pdb_chain(self):
        return self.regions[self.name][0]

    def get_pdb_start(self):
        return self.regions[self.name][1]

    def get_pdb_end(self):
        return self.regions[self.name][2]

    def set_gene(self, gene):
        self.gene = gene

    def add_residue_record(self, name, num):
        self.residue_records[num]=ResidueRecord(name, num)


class FrameworkRegions:
    def __init__(self, chain):
        self.chain = chain
        f1 = [1, 23]
        f1_2 = [43, 56]

        if chain == 'L':
            f2_3 = [73, 106]
        elif chain == 'H':
            f2_3 = [70, 106]
        else:
            sys.exit("Unrecognized chain "+chain)

        f3 = [139, 148]

        self.regions = defaultdict()
        self.regions["F1"] = f1
        self.regions["F1_2"] = f1_2
        self.regions["F2_3"] = f2_3
        self.regions["F3"] = f3
        self.region_names = ["F1", "F1_2","F2_3", "F3"]

    ##Refactor this to use region classes if possible
    def get_regions(self):
        """
        Get the regions as an array of tupples
        """
        regions = []
        regions.append(self.F1())
        regions.append(self.F1_2())
        regions.append(self.F2_3())
        regions.append(self.F3())
        return regions

    def get_residue_record_region(self, region_name):

        res1 = ResidueRecord("-", self.get_start(region_name), self.chain)
        res2 = ResidueRecord("-", self.get_stop(region_name), self.chain)

        region_class = ResidueRegion(res1, res2, region_name)
        return region_class

    def get_residue_record_regions(self):

        res_regions = []
        for name in self.region_names:
            res_region = self.get_residue_record_region(name)
            #print repr(res_region)
            res_regions.append(res_region)

        return res_regions

    def get_start_stop(self, region):
        return self.regions[region]

    def get_start(self, region):
        return self.regions[region][0]

    def get_stop(self, region):
        return self.regions[region][1]

    def F1(self):
        return self.regions["F1"]

    def F1_2(self):
        return self.regions["F1_2"]

    def F2_3(self):
        return self.regions["F2_3"]

    def F3(self):
        return self.regions["F3"]






#############################################################FUTURE2#############################################################    
class molecule:
    def __init__(self, name, atomsormolecule):
        pass
        #Were going to take all atoms one by one and create a molecule.
        #We can take two or more molecules and put them together to create a larger molecule.
        #Examples of a molecule:
            #H and L chain are each molecules.
            #Any ligands are molecules.
            #Molecules can covelently bond to other molecules.
            #Bonding can be a rule based on chemistry in the far flung future.
            

class protein:
    def __init__(self):
        pass

class protein_info:
    """
    The protein has protein molecules.  The PDB has a protein.  Need to write this carefully.
    """
    
    def __init__(self, sequence):
        #Were going to take the sequence.  Find the sequence's family.
        pass
        
    def giveStructure(self, pose):
        #Here, we give the sequence a structure.  We compare it to known structure and find it's family.  We parse it's components into domains using Pfam.
        #We get even more information on those domains. 
        pass
    
    def getInfo(self):
        #Here, we parse Uniprot.  We get all the information we possibly can.
        #We want to KNOW the species.  It's job.  It's fold.  It's required ligand.  Everything we can find on it.
        #Parsable.  Knowable.  
        pass
    
    def getPartners(self):
        #Here, we parse ProtCid. We return all known partners.  We try to determine whether the protein is a monomer or dimer, or something else.
        pass
    
    def breakIntoDomains(self):
        #Here, we break the structure into all domains.  This is useful if we want to do things to the domains.  OR start engineering them.
        pass
    def attachDomains(self, domain1, domain2):
        #Here we attach the two domains.  We will require more information about the design, but this will be added on later.
        pass
    
    
    
#If I ever need it:

class Atom:
    def __init__(self, name):
        pass
    def __add__(self, atom):
        """
        If we add two atoms together.  We create - A molecule, and a bond.
        This would be cool if we could get it to work.
        """
        return Bond(self, atom)
        #self.vdw
        #self.electrons
        #self.energy
        #self.
        pass
    def __iadd__(self):
        return Bond(self, Atom(self.name))
        
        
class Bond:
    def __init__(self, atom1, atom2):
        pass
    
    
    

#Just for Fucks sake:
class nucleus:
    pass
class electron:
    pass

    
