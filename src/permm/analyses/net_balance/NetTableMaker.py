from permm.core.Reaction import Stoic
from warnings import warn

def NetTable(rxn_array, label='Unspecified'):
    result = '%s\n' % (label,)
    result += str(rxn_array)
    return result

def NetTables(mech):
    result = ''
    nrxn_labels = mech.nreaction_dict.keys()
    nrxn_labels.sort()
    for nrxn_label in nrxn_labels:
        nrxn = mech.nreaction_dict[nrxn_label]
        result += NetTable(nrxn, nrxn_label)
        result += '\n'
    
    return result
    