#Jared Adolf-Bryfogle

class SequenceInfo:
    """
    Simple class for holding + accessing sequence metadata

    Original class for sequence info.  Basically deprecated by SequenceStats and PDBConsensusInfo.
    """

    def __init__(self):
        self.start = None
        self.end = None
        self.chain = None

        self.sequence = None
        self.region = None
        self.pdbid = None
        self.pdbpath = None

    def get_sequence(self):
        return self.sequence

    def get_length(self):
        return len(self.sequence)

    def get_pdbID(self):
        return self.pdbID

    def get_pdbpath(self):
        return self.pdbpath

    def get_region(self):
        return self.region

    def get_start_residue(self):
        return self.start

    def get_end_residue(self):
        return self.end

    def get_chain(self):
        return self.chain

    def get_residue(self, resnum):
        """
        If region is given, resnum is residue number of PDB
        If not, resnum in Rosetta resnum
        """
        print self.sequence
        print repr(resnum)
        if self.start:
            index_num = resnum-self.start
            one_letter_code = self.sequence[index_num]

        else:
            one_letter_code = self.sequence[resnum-1]

        return one_letter_code

    def set_sequence(self, sequence):
        self.sequence = sequence
    def set_pdbID(self, pdbID):
        self.pdbID = pdbID
    def set_pdbpath(self, pdbpath):
        self.pdbpath = pdbpath
    def set_region(self, region):
        self.region = region
        rSP = region.split(":")
        self.start = int(rSP[0])
        self.end = int(rSP[1])
        self.chain = rSP[2]
