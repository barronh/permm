__all__ = ['IRRTable']

def IRRTable(mech, delim = '_'):
    reaction_keys = [(int(k.split(delim)[1]),k) for k in mech.reaction_dict.keys()]
    reaction_keys.sort()
    reaction_keys = [k[1] for k in reaction_keys]
    result = ''
    template = '%s\t%13.7f\t%s\n'
    for rkey in reaction_keys:
        result += template % (rkey, mech.irr[rkey].sum(), mech.nreaction_dict[rkey].sum())
        
    return result