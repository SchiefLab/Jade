


class SequenceResults:
    """
    Simple class for holding, calculating, + accessing result data
    Residue Numbers are in Rosetta numbering.

    Original class for sequence stats.  Basically deprecated by SequenceStats and PDBConsensusInfo.
    """
    def __init__(self):
        self.data = dict()
        self.reference = dict()

    def add_residue(self, resnum, one_letter_code, decoy):
        if not self.data.has_key(resnum):
            self.data[resnum]=dict()
            self.data[resnum][one_letter_code]=dict()
            self.data[resnum][one_letter_code]['freq']=1
            self.data[resnum][one_letter_code]['decoys']=[]
            self.data[resnum][one_letter_code]['decoys'].append(decoy);#This is to keep track of which decoys have which mutation.

        else:
            if not self.data[resnum].has_key(one_letter_code):
                self.data[resnum][one_letter_code]=dict()
                self.data[resnum][one_letter_code]['freq']=0
                self.data[resnum][one_letter_code]['decoys']=[]
            self.data[resnum][one_letter_code]['freq']+=1
            self.data[resnum][one_letter_code]['decoys'].append(decoy)

    def add_reference_residue(self, resnum, one_letter_code):
        self.reference[resnum]=one_letter_code

    def get_freq(self, resnum, one_letter_code):
        try:
            x = self.data[resnum][one_letter_code]['freq']
            return x
        except KeyError:
            return 0

    def get_total(self, resnum):
        total = 0
        for code in self.data[resnum]:
            freq = self.get_freq(resnum, code)
            total = total +freq
        return total

    def get_percent(self, resnum, one_letter_code):
        total = self.get_total(resnum)
        freq = self.get_freq(resnum, one_letter_code)

        percent = float(freq)/float(total)
        return percent

    def get_percent_string(self, resnum, one_letter_code):
        return "%.2f"%self.get_percent(resnum, one_letter_code)

    def get_reference_residue(self, resnum):
        return self.reference[resnum]

    def get_all_residues_observed(self, resnum):
        return sorted(self.data[resnum].keys())

    def get_all_residue_numbers(self):
        return sorted(self.data.keys())

    def get_decoys_with_aa(self, resnum, one_letter_code):
        """
        Returns all decoys with a specific mutation at a position.
        """
        try:
            return self.data[resnum][one_letter_code]['decoys']
        except KeyError:
            return []

    def get_decoys_with_joint_aa(self, resnum_one_letter_code_pair):
        """
        Will output decoys that have x, y, z mutations at positions a, b, c
        """
        pass

    ### reference Comparison Functions ###
    def get_all_mutated_positions(self):
        mutated_positions = []
        for resnum in self.data:
            if not self.reference.has_key(resnum):
                print "Position in data does not match position in reference"

            if self.get_percent(resnum, self.reference[resnum])==1.0:
                pass
            else:
                mutated_positions.append(resnum)

        if mutated_positions:return mutated_positions
        else:print "No mutations found"

    def get_all_reference_percent_observed(self):
        """
        Returns array of tripplets of [postion, one_letter_code, percent] of reference amino acid found.
        """
        tripplet_array = []
        for resnum in self.reference:
            if not self.data.has_key(resnum):
                print "Position in reference does not match any position in data"

            percent = self.get_percent(resnum, self.reference[resnum])
            tripplet = [resnum, self.reference[resnum], percent]
            tripplet_array.append(tripplet)
        return tripplet_array
