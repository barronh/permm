__all__ = ['PhyTable']
from numpy import where
import warnings
from permm.core.Species import Species
warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: '%s: %s\n' % ('ProcessWarning', message)


default_procs = ['Initial', 'Emissions', 'Chemistry', 'H_Trans', 'V_Trans', 'Entrain', 'Deposit', 'TEMPADJ', 'Final', 'Proc_Sum']
default_proc_split = ['Chemistry', 'H_Trans', 'V_Trans', 'Entrain']

def PhyTables(mech):
    result = ''
    for spck, spcv in mech.species_dict.iteritems():
        if isinstance(spcv, Species):
            result += PhyTable(mech, spck)
            result += '\n'
    
    return result
    
def VOCTable(mech):
    result = PhyTable(mech, 'VOC', speciated = True)
    result += PhyTable(mech, 'VOCC', speciated = True)
    return result
    
def PhyTable(mech, spc, extended = True, process_list = default_procs, split_list = default_proc_split, speciated = False, initial = 'Initial'):
    tag = ('', '  --XT')[extended]
    
    result = '%s%s\n' % (spc, tag)
    times = tuple(mech.mrg.variables['TFLAG'][:,0,1].tolist())
    ntimes = len(times)
    result += ('  %-19s'+"%11.4i"*ntimes+'%11s\n') % (('Time',)+times+('Sum',))
    
    if speciated:
        spc_grp = eval(spc, None, mech.species_dict)
        spcs = spc_grp.names()+[spc]
    else:
        spcs = [spc]
        
    for process in process_list:
        if speciated:
            result += '  %s\n' % process
            
        for spc in spcs:
            if speciated:
                label = "  %s" % (spc,)
            else:
                label = process
            
            if process == initial:
                result += display(mech, spc, process, label, lambda x: x, lambda x: x[0])
            else:
                result += display(mech, spc, process, label, lambda x: x)

        if extended and process in split_list:
            if speciated:
                result += '   Gains\n'

            for spc in spcs:
                if speciated:
                    label = "  %s" % (spc,)
                else:
                    if extended:
                        label = '%s, gain' % (process,)
                    else:
                        label = process
                        
                result += display(mech, spc, process, label, lambda x: where(x > 0, x, 0))

            if speciated:
                result += '   Losses\n'
            
            for spc in spcs:
                if speciated:
                    label = "  %s" % (spc,)
                else:
                    if extended:
                        label = '%s, loss' % (process,)
                    else:
                        label = process
                        
                result += display(mech, spc, process, label, lambda x: where(x < 0, x, 0))

    return result
    
def display(mech, spc, process, label, condition = None, agg = lambda x: x.sum(), chemdetail = False):
    if condition is None:
        condition = lambda x: x
    try:
        hourly = condition(mech("%s[%s].array()" % (process, spc)))
        daily = hourly.sum()
        n_vals = hourly.size+1
        values = hourly

        display_values = (label,)+tuple(values.tolist())+(agg(values),)
        result = ("  %-19s"+"%11.4f"*n_vals+'\n') % display_values
    except KeyError, strerror:
        result = ''
        warnings.warn('%s' % (strerror))

    return result