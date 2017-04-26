from collections import defaultdict

class CycleData:
    """
    Simple class to hold data.  Should be replaced.
    """
    def __init__(self, pdb_path, pdb_id):
        self.pdb_path = pdb_path
        self.pdb_id = pdb_id
        self.values = defaultdict(dict)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, item):
        return self.get_cycle(int(item))

    def __repr__(self):
        return repr(self.values)

    def get_cycle(self, n):
        return self.values[int(n)]


    def get_acceptance(self, n):
        return self[n]["mc_accept"]

    def get_mc_e(self, n):
        return self[n]["mc_e"]

    def get_final_e(self,n):
        return self[n]["final_e"]

    def get_protocol_final_e(self):
        return self.final_e

    def get_native_e(self):
        return self.get_protocol_final_e()

    def set_protocol_final_e(self, e):
        self.final_e = float(e)

    def set_protocol_start_e(self, e):
        self.start_e = float(e)

    def set_protocol_native_e(self, e):
        self.native_e = float(e)


    def set_final_energy(self, n, e):
        n = int(n)
        e = float(e)
        self.values[n]["final_e"] = e

    def set_mc_data(self, n, e, acceptance):
        n = int(n)
        e = float(e)
        #acceptance = bool(acceptance)
        self.values[n]["mc_e"] = e
        self.values[n]["mc_accept"] = int(acceptance)

    def set_cum_acceptance(self, n, cum_accepts):
        self.values[n]["cum_accepts"] = cum_accepts
        print "Setting "+repr(n)+" "+repr(cum_accepts)

    def get_cum_acceptance(self, n):
        return self.values[n]["cum_accepts"]

    def set_delta_e(self, n, delta_e):
        self.values[n]["delta_e"] = delta_e

    def get_delta_e(self, n):
        return self.values[n]["delta_e"]