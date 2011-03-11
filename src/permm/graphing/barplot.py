__all__ = ['profile_process']

from warnings import warn
from pylab import figure, title, bar, savefig, legend, grid, xlabel, ylabel, xticks, axis
from matplotlib.cm import get_cmap
from matplotlib.colors import rgb2hex
from numpy import arange, ndarray

def profile_process(mech, species_group, process, loc = 1, width = .8, sort = False, fig = None, cmap = None, normalize = True, missing = None):

    # If necessary, start a new figure.  Otherwise,
    # store the previous colors as a dictionary
    if fig is None:
        fig = figure()
        used_colors = {}
    else:
        used_colors = dict([(p.get_label(),p.get_facecolor()) for p in fig.axes[0].patches])
    # Evaluate species group name
    # and return a species object
    species_group = mech(species_group)

    # Evaluate process group name
    # and return a process object
    process_obj = mech(process)
    
    # Get the list of component species
    component_species = species_group.names()
        
    # Bar chart starts at the bottom
    bottom, top = 0, 0
    
    # Create a receptacle for values
    vals = []
    
    # For each component species, get the process value
    # Use zero where missing
    for species in component_species:
        factor = species_group[species]
        species_obj = mech(species)
        try:
            vals.append((species, float(factor) * process_obj[species_obj].array().sum()))
        except KeyError, (m):
            warn(m)
            if missing is not None:
                vals.append((species, missing))

    # The number of patches is based on the number of components
    nbars = len(vals)
    colors = iter(get_cmap(cmap)(arange(nbars, dtype = 'f')/(nbars-1)))

    # Normalize to 100% if requested
    if normalize:
        total = process_obj[species_group].array().sum()
        vals = [ (k, v/total) for k, v in vals]

    # Sort based on value if requested
    if sort:
        vals = [ (v, k) for k, v in vals]
        vals.sort(reverse = True)
        vals = [ (k, v) for v, k in vals]

    # Plot each component value using pre-existing
    # color where possible
    for name, val in vals:
        if used_colors.has_key(name):
            color = used_colors[name]
        else:
            while True:
                color = rgb2hex(colors.next()[:-1])
                if color not in used_colors.values():
                    break
            used_colors[name] = color    
        
        bar(loc, val, width, bottom = bottom, label = name, color = color )
        bottom += val
    return fig
    
if __name__ == '__main__':
    from pyPA.utils.util import AttrDict
    from net_balance import get_pure_mech
    from permm.netcdf import NetCDFFile
    conf = AttrDict()
    conf['mechanism'] = 'cb05'
    conf['model'] = 'camx'
    conf['mrgfile'] = '/Users/barronh/Development/net_balance/src/net_balance/testdata/test.mrg.nc'
    conf['species']  = dict(HNO3=dict(ymin=None),
                        NO2=dict(ymin=None),
                        NO=dict(ymin=None),
                        NOx=dict(ymin=None),
                        NOyN=dict(ymin=None),
                        NTR=dict(ymin=None)
                        )
    conf['process'] = {'H_Trans': dict(color = 'b'), 
                    'V_Trans': dict(color = 'g'),
                    'Emissions': dict(color = 'r'),
                    'Deposit': dict(color = 'c'),
                    'Aero_Chem': dict(color = 'm'),
                    'CHEM': dict(color = 'y'),
                    'Motion': dict(color = 'teal'),
                    'TEMPADJ': dict(color = 'aqua')
                    }
    conf['init'] = 'INIT'
    conf['final'] = 'FCONC'
    mech = get_pure_mech('_'.join([conf['mechanism'],conf['model']]))
    mrg_file = NetCDFFile(conf['mrgfile'],'r')
    mech.set_mrg(mrg_file)

    conf['chem'] = 'CHEM'
    conf['title'] = '%(species)s IPR plot'
    conf['outdir'] = '.'
    conf['end_date'] = True
    fig = profile_process(mech, 'NOyN', 'INIT', loc = 0, width = .2, sort = False, missing = None)
    fig = profile_process(mech, 'NOyN', 'Emissions', loc = 0.2, width = .2, sort = False, fig = fig, missing = None)
    fig = profile_process(mech, 'NOyN', 'H_Trans', loc = 0.4, width = .2, sort = False, fig = fig, missing = None)
    fig = profile_process(mech, 'NOyN', 'V_Trans', loc = 0.6, width = .2, sort = False, fig = fig, missing = None)
    fig = profile_process(mech, 'NOyN', 'FCONC', loc = .8, width = .2, sort = False, fig = fig, missing = None)
    patches={}
    for patch in fig.axes[0].patches:
        patches.setdefault(patch.get_label(), patch)
    names = patches.keys()
    patches = [patches[k] for k in names]
    
    fig.legend(patches, names)
    xticks(arange(5, dtype = 'f')/5+.1, ['Initial', 'Emis', 'Horz', 'Vert', 'Final'])
    axis('tight')
    fig.savefig('NOyN_bar.pdf', fmt = 'pdf')