__all__ = ['irr_plot', 'phy_plot', 'phy_plots']


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
from numpy import array, concatenate, zeros, arange, ceil, concatenate
from pylab import figure, title as pylabtitle, plot_date, savefig, legend, axis, grid, axes, xlabel, ylabel, subplot, gca, twinx, rcParams
from matplotlib.dates import DateFormatter, date2num
from matplotlib.cm import get_cmap
from matplotlib.font_manager import FontProperties
import re
import operator
import os

if rcParams['text.usetex']:
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

def get_dates(mrg_file, end_date = True):
    date_ints = mrg_file.variables['TFLAG'][...,0,:].reshape(-1, 2)
    date_objs = array([datetime.strptime("%iT%06i" % (d,t), "%Y%jT%H%M%S") for d,t in date_ints]).reshape(mrg_file.variables['TFLAG'].shape[:-2])
    if end_date:
        date_objs = concatenate([date_objs[[0]]-(date_objs[-1]-date_objs[-2]), date_objs])
    else:
        date_objs = concatenate([date_objs, date_objs[[-1]]+(date_objs[-1]-date_objs[-2])])

    return date_objs
    
def get_date_steps(mrg_file, end_date = True):
    date_objs = get_dates(mrg_file, end_date)

    date_objs = date_objs.repeat(2,0)[1:-1]
    return date_objs

def none_defaults_to(v, default):
    return {True: default, False: v}[v is None]

def plot(mech, y, stepped = True, end_date = True, figure_settings = {}, axis_settings = {}, line_settings = {}):
    
    if stepped:
        date_objs = get_date_steps(mech.mrg, end_date)
        y = y.repeat(2,0)
    else:
        date_objs = get_dates(mech.mrg, end_date)

    axis_settings = axis_settings.copy()
    axis_settings.setdefault('xlabel', 'Time')
    line_settings = line_settings.copy()
    line_settings.setdefault('linestyle', '-')
    line_settings.setdefault('linewidth', 3)
    line_settings.setdefault('marker', 'None')

    fig = figure(**figure_settings)
    ax = fig.add_subplot(111, **axis_settings)
    
    ax.plot_date(date_objs, y, **line_settings)
    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    fig.autofmt_xdate()

    return fig
    
def irr_plot(
             mech, reactions, species,
             factor = 1, units = None, end_date = True,
             chem = 'CHEM', title = None, ylim = None, xlim = None, 
             figure_settings = {}, axis_settings = {}, line_settings = {},
             fig = None, cmap = None, ncol = 2
            ):
    """
    irr_plot makes a plot with lines for each reaction object provided
    
      mech - Mechanism object
      reactions - list of reaction objects
      species - species to plot from reactions
      end_date - dates are provided for the interval end
      units - string for y-axis label (defaults to mech.irr.units)
      ncol - number of columns in the legend
      fig - a previous figure to plot results on
      cmap - color map to be used
      chem - the process name of the total chemistry or None to 
             not plot total chemistry
      factor - a scaling factor for process data (i.e. if data is ppm 
               and you want to plot ppt, factor = 1e6)
      title - a string for the title of the figure
      ylim - an iterable (e.g. tuple or list) of min and max values
             for y-scale
      xlim - an iterable (e.g. tuple or list) of min and max values
             for x-scale
      figure_settings - properties for the figure; dictionary of 
                        keywords that pertain to the matplotlib
                        figure function
      axis_settings - properties for the figure axes; dictionary of 
                      keywords that pertain to the matplotlib
                      axes function
      line_settings - properties for all lines; dictionary of keywords
                      for matplotlib plot_date function
    """
    
    date_objs = get_date_steps(mech.mrg, end_date)
    units = none_defaults_to(units, mech.irr.units)
    nlines = len(reactions)
    if rcParams['text.usetex']:
        title_str = 'Plot of \ce{%s} Reactions' % (species.name, )
    else:
        title_str = 'Plot of %s Reactions' % (species.name, )

    title_str = none_defaults_to(title, title_str)
    
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    if fig is None:
        fig = figure(**figure_settings)
        legend_line_height = .037
        if rcParams['text.usetex']: legend_line_height = .045
        nrows = ceil(nlines/ncol)
        legend_space = (nrows + 2) * legend_line_height
        ax = axes([.1, .1 + legend_space, .8, .8 - legend_space], **axis_settings)
    else:
        ax = gca()
        
    grid(True)
    pylabtitle(title_str % locals())

    options = line_settings
    options.setdefault('linestyle', '-')
    options.setdefault('linewidth', 3)
    options.setdefault('marker', 'None')
    for rxn in reactions:
        options['color'] = colors.next()
        data = rxn[species] * factor
        rxn_sum = rxn.sum()
        reaction_label = rxn.display(digits = None)
        options['label'] = reaction_label
        if rcParams['text.usetex']: options['label'] = '\ce{' + options['label'].replace('>', ']').replace('=', '->[') + '}'
        
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
    if ylim is not None:
        ax.set_ylim(*ylim)

    if xlim is not None:
        ax.set_xlim(*xlim)

    fig.autofmt_xdate()
    labels = [l.get_label() for l in ax.lines]
    handles = ax.lines
    leg = fig.legend(handles, labels, ncol = ncol, loc = 8, prop = FontProperties(size=10))

    return fig

def phy_plot(mech, species, init = 'INIT', final = 'FCONC', factor = 1, end_date = True,
             exclude = ['UCNV', 'AVOL', 'DATE', 'TIME', 'I', 'J', 'K'], processes = None,
             units = None, ylim = None, xlim = None,
             title = None, filter = True, fig = None, ncol = 1, cmap = None,
             figure_settings = {}, axis_settings = {},
             line_settings = {'linestyle': '-', 'linewidth': 3, 'marker': 'None'},
             **kwds):

    """
    mech - perm.Mechanism.Mechanism object
    species - perm.SpeciesGroup.Species object
    title - title
    init - Name of initial concentration process
    final - Name of final concentration process
    ncol - number of legend columns (default: 1)
    fig - figure to plot on (default: None)
    cmap - matplotlib color map for lines (default: None)
    filter - remove processes with zero values (default: True)
    kwds - * <process name1> - process names from mech.process_dict can be 
                               provided to limit the processes shown.  When 
                               provided, a process should be a dictionary of 
                               matplotlib plot options (common: linestyle,
                               linewidth, label, marker).  The dictionary can
                               be empty.
           * <process name2> - same as process 1
           * <process nameN> - same as process 1
           * CONC - same as process 1, but for the concentration line
           * end_date - times are for the time period end
           
           
    """
    processes = none_defaults_to(processes, set(kwds.keys()).intersection(mech.process_dict.keys()))
    if processes == set():
        processes = set([k for k, v in mech.process_dict.iteritems() if len(v.names)==1 and v.names[0]==k]).difference([init, final]+exclude)

    processes = dict(zip(processes, map(lambda k: kwds.get(k,line_settings.copy()), processes)))
            
    if rcParams['text.usetex']:
        title_str = '\ce{%s} Processes' % species.name
    else:
        title_str = '%s Processes' % species.name

    title_str = none_defaults_to(title, title_str)

    if isinstance(processes, (list, set)):
        processes = dict([(k,{}) for k in processes])

    units = none_defaults_to(units, mech.ipr.units)
    nlines = len(processes)
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    
    date_objs = get_date_steps(mech.mrg, end_date)
    legend_width = ncol * max([len(k) for k in processes.keys()]) * 10 # 10 is legend.fontsize
    
    if fig is None:
        fig = figure(figsize = (8,6))
        legend_pct = legend_width / fig.bbox.get_points()[1,0]
        ax = axes([0.1, 0.12, 0.8 - legend_pct, 0.8])
    else:
        ax = gca()
    
    tax = ax
    grid(True)
    pylabtitle(title_str % locals())
    options = kwds.get('CONC',line_settings.copy())
    options.setdefault('color', 'k')
    options.setdefault('label', 'Conc')
    options.setdefault('linestyle', '-')
    options.setdefault('linewidth', 3)
    options.setdefault('marker', 'None')
    init_vals = mech('%s'  % (init,))[species].array()
    fconc_vals = mech('%s'  % (final,))[species].array()
    data = concatenate([init_vals[..., :1], fconc_vals[...]], axis = -1) * factor
    tax.plot_date(get_dates(mech.mrg), data, **options)
    
    for process in processes.keys():
        options = processes[process]
        options.setdefault('color', colors.next())
        options.setdefault('linestyle', '-')
        options.setdefault('linewidth', 3)
        options.setdefault('marker', 'None')
        process_label = process
        if rcParams['text.usetex']:
            process_label = process_label.replace('_', '\_')
            
        options.setdefault('label', process_label)
        data = mech('(%s)' % (process,))[species].array().repeat(2,0) * factor
        if data.nonzero()[0].any() or not filter:
            plot_date(date_objs, data, **options)
            
    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    ax.set_xlim(date2num(date_objs[0]), date2num(date_objs[-1]))
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)

    fig.autofmt_xdate()
    legend(ncol = ncol, loc = (1.05,0), prop = FontProperties(size=10))
    return fig

def phy_plots(conf, filter = True, fmt = 'pdf', fig_append = 'IPR'):
    add_mech(conf)
    
    mech = conf['mech']
    date_objs = get_date_steps(mech.mrg, conf)

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
