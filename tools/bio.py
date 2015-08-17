
## Biopython functions

def has_id(model, id):
    """
    Returns true or false if the model has the chain.  Because biopython is not updating it's index that has_id is using.  WTF.
    """
    for i in model:
        if i.id == id:
            return True
    return False


def get_chain_length(bio_chain):

    l = 0
    for res in bio_chain:
        if res.id[0]==' ':
            l+=1
    return l

def get_biochain_sequence(bio_chain):
    seq = ""
    d = definitions()

    for res in bio_chain:
        if res.id[0]==' ':
            aa = d.get_one_letter_from_three(res.resname)
            if not aa:
                print "Skipping non-canonical resname: "+res.resname
                print "This could pose a problem!"
                continue
            seq = seq+aa
    return seq
