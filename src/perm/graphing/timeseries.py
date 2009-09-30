__all__ = ['irr_plot', 'irr_plots', \
           'loss_plots', 'phy_plots', \
           'prod_plots', 'rxn_plots']


HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from warnings import warn
from datetime import datetime
from types import InstanceType
from ..netcdf import NetCDFFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from numpy import array, concatenate, zeros, arange, ceil
from pylab import figure, title, plot_date, savefig, legend, axis, grid, axes, xlabel, ylabel, subplot, gca, twinx, rcParams
from matplotlib.dates import DateFormatter, date2num
from matplotlib.cm import get_cmap
from matplotlib.font_manager import FontProperties
import re
import operator
import os

if rcParams['text.usetex'] == True:
    rcParams['text.latex.preamble'] = '\usepackage[font=sf]{mhchem}'

def add_mech(conf):
    if conf.has_key('mech'):
        mech = conf['mech']
    else:
        from net_balance import get_pure_mech
        mech = conf['mech'] = get_pure_mech('_'.join([conf['mechanism'].lower(), conf['model'].lower()]))
        if isinstance(conf['mrgfile'], (PseudoNetCDFFile, InstanceType)):
            mrg_file = conf['mrgfile']
        else:
            mrg_file = NetCDFFile(conf['mrgfile'],'r')

        mech.set_mrg(mrg_file)

def get_date(mrg_file, conf):
    date_ints = mrg_file.variables['TFLAG'][:,0,:]
    date_objs = array([datetime.strptime("%iT%06i" % (d,t), "%Y%jT%H%M%S") for d,t in date_ints])
    if conf.get('end_date',True):
        date_objs = concatenate([date_objs[[0]]-(date_objs[-1]-date_objs[-2]), date_objs])
    else:
        date_objs = concatenate([date_objs, date_objs[[-1]]+(date_objs[-1]-date_objs[-2])])

    date_objs = date_objs.repeat(2,0)[1:-1]
    return date_objs

def irr_plot(mech, reactions, species, **conf):
    """
    irr_plot makes a plot with lines for each reaction
      conf - configuration dictionary
      mech - Mechanism object
      date_objs - dates assoicated with IRR
      species - species name for production plot
      nlines - limit number of process lines to nlines
      combine - iterable of iterables; inner iterable should be used to make net reactions
      fig - figure to plot on
      cmap - color map to be used
    """
    combine = conf.get('combine',[()])
    reactions = [ rxn for rxn in reactions if rxn not in reduce(operator.add, combine) ]
    date_objs = get_date(mech.mrg, conf)
    units = conf.get('units', mech.irr.units)
    nlines = min(conf.get('nlines',8),len(reactions)+1)
    ncol = float(conf.get('ncol',2))
    fig = conf.get('fig', None)
    cmap = conf.get('cmap', None)
    if rcParams['text.usetex']:
        title_str = 'Plot of \ce{%s} for %d Reactions' % (species.name, len(reactions))
    else:
        title_str = 'Plot of %s for %d Reactions' % (species.name, len(reactions))

    title_str = conf.get('title', title_str)
    chem = conf.get('chem', 'Chemistry')
    factor = conf.get('factor', 1.)
    
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    if fig is None:
        fig = figure()
        ax = axes([.1, .1+.037*(ceil(nlines/ncol)+1), .8, .8-.037*ceil(nlines/ncol)])
    else:
        ax = gca()
        
    grid(True)
    title(title_str % locals())

    if combine != [()]:
        reactions = reactions + map(lambda t2: '+'.join(t2), combine)

    reactions = [ (abs(mech('(%s)' % (rxn))[species]).sum(),rxn) for rxn in reactions]

    reactions.sort(reverse = True)

    reactions = [r for v,r in reactions]
    
    options = {}
    options.setdefault('linestyle', '-')
    options.setdefault('linewidth', 3)
    options.setdefault('marker', 'None')

    for rxn in reactions[:nlines-1]:
        options['color'] = colors.next()
        data = mech('(%s)' % (rxn, ))[species] * factor
        options['label'] = re.compile('\d+.\d+\*').sub('', str(mech(rxn).sum())).replace(' ', '')
        plot_date(date_objs, data.repeat(2,0), **options)

    data = zeros(data.shape, dtype = data.dtype)
    for rxn in reactions[nlines-1:]:
        data += mech('(%s)' % (rxn,))[species]
        
    if (data != 0).any():
        options['label'] = 'other'
        plot_date(date_objs, data.repeat(2,0), **options)
    if chem is not None:
        options['color'] = 'black'
        options['marker'] = 'o'
        options['label'] = 'Total Chem'
        try:
            data = mech('%s'  % (chem,))[species].array()
        except:
            warn('Using sum of reactions for %(species)s' % locals())
            data = mech.make_net_rxn(species, species, False)[species]
    
        plot_date(date_objs, data.repeat(2,0), **options)

    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    ax.set_xlim(date2num(date_objs[0]), date2num(date_objs[-1]))
    if conf.has_key('ylim'):
        ax.set_ylim(*conf['ylim'])

    fig.autofmt_xdate()
    labels = [l.get_label() for l in ax.lines]
    handles = ax.lines
    leg = fig.legend(handles, labels, ncol = ncol, loc = 8, prop = FontProperties(size=10))

    return fig

def phy_plot(mech, species, **kwds):
    """
    mech - perm.Mechanism.Mechanism object
    species - perm.SpeciesGroup.Species object
    kwds - * title - title
           * init - Name of initial concentration process
           * final - Name of final concentration process
           * linestyle - default line style (default: '-')
           * linewidth - default line width (default: 3)
           * marker - default line marker style (default: None)
           * ncol - number of legend columns (default: 1)
           * fig - figure to plot on (default: None)
           * cmap - matplotlib color map for lines (default: None)
           * filter - remove processes with zero values (default: True)
           * <process name1> - process names from mech.process_dict can be 
                               provided to limit the processes shown.  When 
                               provided, a process should be a dictionary of 
                               matplotlib plot options (common: linestyle,
                               linewidth, label, marker).  The dictionary can
                               be empty.
           * <process name2> - same as process 1
           * <process nameN> - same as process 1
           * end_date - times are for the time period end
           
           
    """
    init = kwds.get('init', 'INIT')
    final = kwds.get('final', 'FCONC')
    factor = kwds.get('factor', 1)
    exclude = kwds.get('exclude', 'UCNV AVOL DATE TIME I J K'.split())
    if kwds.has_key('processes'):
        processes = kwds['processes']
    else:
        processes = set(kwds.keys()).intersection(mech.process_dict.keys())
        if processes == set():
            processes = set([k for k, v in mech.process_dict.iteritems() if len(v.names)==1 and v.names[0]==k]).difference([init, final]+exclude)

    processes = dict(zip(processes, map(lambda k: kwds.get(k,{}), processes)))
            
    filter = kwds.get('filter', True)
    fig = kwds.get('fig', None)
    ncol = kwds.get('ncol', 1)
    if rcParams['text.usetex']:
        title_str = '\ce{%s} Processes' % species.name
    else:
        title_str = '%s Processes' % species.name
    title_str = kwds.get('title', title_str)
    cmap = kwds.get('cmap', None)
    linestyle = kwds.get('linestyle', '-')
    linewidth = kwds.get('linewidth', 3)
    marker = kwds.get('marker', 'None')

    if isinstance(processes, (list, set)):
        processes = dict([(k,{}) for k in processes])

    units = kwds.get('units', mech.ipr.units)
    nlines = len(processes)
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    
    date_objs = get_date(mech.mrg, kwds)
    legend_width = ncol * max([len(k) for k in processes.keys()]) * 10 # 10 is legend.fontsize
    
    if fig is None:
        fig = figure(figsize = (8,6))
        legend_pct = legend_width / fig.bbox.get_points()[1,0]
        ax = axes([0.1, 0.12, 0.8 - legend_pct, 0.8])
    else:
        ax = gca()
    
    tax = ax
    grid(True)
    title(title_str % locals())
    options = {'color': 'k'}
    options.setdefault('linestyle', linestyle)
    options.setdefault('linewidth', linewidth)
    options.setdefault('marker', 'o')
    data = mech('%s'  % (init,))[species].array() * factor
    tax.plot_date(date_objs[::2], data, **options)
    options.setdefault('label', 'Conc')
    options['marker'] = 'x'
    data = mech('%s'  % (final,))[species].array() * factor
    tax.plot_date(date_objs[1::2], data, **options)
    
    for process in processes.keys():
        options = processes.get(process, {})
        options.setdefault('color', colors.next())
        options.setdefault('linestyle', linestyle)
        options.setdefault('linewidth', linewidth)
        options.setdefault('label', process)
        options.setdefault('marker', marker)
        data = mech('(%s)' % (process,))[species].array().repeat(2,0) * factor
        if data.nonzero()[0].any() or not filter:
            plot_date(date_objs, data, **options)
            
    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    ax.set_xlim(date2num(date_objs[0]), date2num(date_objs[-1]))
    if kwds.has_key('ylim'):
        ax.set_ylim(*kwds['ylim'])

    fig.autofmt_xdate()
    legend(ncol = ncol, loc = (1.05,0), prop = FontProperties(size=10))
    return fig

def phy_plots(conf, filter = True, fmt = 'pdf', fig_append = 'IPR'):
    add_mech(conf)
    
    mech = conf['mech']
    date_objs = get_date(mech.mrg, conf)

    for species, species_options in conf['species'].iteritems():
        try:
            fig = phy_plot(conf, mech,  date_objs, species, species_options, filter = filter)
            fig.savefig(os.path.join(conf['outdir'], '%s_%s.%s' % (species, fig_append, fmt)), format = fmt)
        except KeyError, detail:
            warn(detail)
    
if __name__ == '__main__':
    conf = dict()
    conf['mechanism'] = 'cb05'
    conf['model'] = 'camx'
    conf['mrgfile'] = '/Users/barronh/Development/net_balance/src/net_balance/testdata/test.mrg.nc'
    conf['species']  = dict(HNO3={},
                        NO2={},
                        NO={},
                        NOx={},
                        NOyN={},
                        NTR={}
                        )
    conf['species'] = dict(NOyN={})
    conf['process'] = {'H_Trans': {}, 
                    'V_Trans': {},
                    'Emissions': {},
                    'Deposit': {},
                    'Aero_Chem': {},
                    'CHEM': {},
                    'Motion': {},
                    'TEMPADJ': {}
                    }
    conf['init'] = 'INIT'
    conf['final'] = 'FCONC'
    conf['chem'] = 'CHEM'
    conf['title'] = '%(species)s IPR plot'
    conf['outdir'] = '.'
    conf['end_date'] = True
    phy_plots(conf)
#    rxn_plot(conf, combine = [('RXN_01', 'RXN_02', 'RXN_03')])
    conf['title'] = '%(species)s IRR plot'
    rxn_plots(conf)
    conf['title'] = '%(species)s Chem plot'
    chem_plots(conf)
