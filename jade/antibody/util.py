import sys

from jade.basic.structure import Structure


def get_overhang_sele(pymol_writer, cdr, overhang =3):
    """
    Gets the selection for terminal superposition in PyMol
    """
    if not isinstance(cdr, Structure.CDR): sys.exit()
    start = cdr.get_pdb_start()
    end = cdr.get_pdb_end()

    chain = cdr.get_pdb_chain()

    start_sele = start - overhang
    end_sele = end + overhang

    resi = []
    resi.append((start_sele, start))
    resi.append((end, end_sele))

    sele = pymol_writer.get_sele(chain, resi)
    print sele
    return sele

def get_all_cdr_sele(cdr, stem = 0):
    """
    If any stem residues are given, will include stem residues in selection.
    """
    if not isinstance(cdr, Structure.CDR): sys.exit()

    #chain L and resid 24-42
    sele = " ".join(["chain", cdr.get_pdb_chain(), "and", "resid", str((cdr.get_pdb_start() - stem)) + "-" + str((cdr.get_pdb_end() + int(stem))) ])
    print sele
    return sele
