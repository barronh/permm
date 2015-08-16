__all__ = ['irr_plot', 'phy_plot', 'phy_plots']


HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from warnings import warn
from datetime import datetime, timedelta
from types import InstanceType
from permm.netcdf import NetCDFFile
from PseudoNetCDF.sci_var import PseudoNetCDFFile
from numpy import array, concatenate, zeros, arange, ceil, concatenate, diff
from pylab import figure, title as pylabtitle, plot_date, savefig, legend, axis, twinx, rcParams
from matplotlib.dates import DateFormatter, date2num
from matplotlib.cm import get_cmap
from matplotlib.font_manager import FontProperties
import re
import operator
import os

if rcParams['text.usetex']:
    rcParams['text.latex.preamble'] = '\usepackage[font=sf]{mhchem}'

def get_dates(mrg_file, end_date = True, slice = slice(None)):
    if 'TFLAG' in mrg_file.variables:
        date_ints = mrg_file.variables['TFLAG'][...,0,:].reshape(-1, 2)
        date_objs = array([datetime.strptime("%iT%06i" % (d,t), "%Y%jT%H%M%S") for d,t in date_ints]).reshape(mrg_file.variables['TFLAG'].shape[:-2])[slice]
    elif 'time' in mrg_file.variables or 'tau0' in mrg_file.variables:
        try:
            date_var = mrg_file.variables['time']
            unit, start = map(lambda x: x.strip(), date_var.units.split('since'))
            end_date = False
        except:
            date_var = mrg_file.variables['tau0']
            end_date = False
            unit, start = map(lambda x: x.strip(), date_var.units.split('since'))
        try:
            sdate = datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
        except:
            sdate = datetime.strptime(start, "%Y-%m-%d %H:%M:%S %Z")
        date_objs = array([sdate + timedelta(**{unit: i}) for i in date_var[:].tolist()])
    if end_date:
        date_objs = concatenate([date_objs[[0]]-(date_objs[-1]-date_objs[-2]), date_objs])
    else:
        date_objs = concatenate([date_objs, date_objs[[-1]]+(date_objs[-1]-date_objs[-2])])

    return date_objs
    
def get_date_steps(mrg_file, end_date = True, slice = slice(None)):
    date_objs = get_dates(mrg_file, end_date, slice = slice)

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
             factor = 1, slice = slice(None), units = None, end_date = True,
             chem = 'CHEM', title = None, ylim = None, xlim = None, 
             figure_settings = {}, axis_settings = {}, line_settings = {},
             fig = None, cmap = None, ncol = 2
            ):
    __doc__ = """
    irr_plot makes a plot with lines for each reaction object provided
    
      mech - Mechanism object
      reactions - list of reaction objects
      species - species to plot from reactions
      end_date - dates are provided for the interval end
      units - string for y-axis label (defaults to mech.irr_dict.units)
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
    
    date_objs = get_date_steps(mech.mrg, end_date, slice = slice)
    units = none_defaults_to(units, getattr(mech.irr_dict, 'units', 'Unknown'))
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
        ax = fig.add_axes([.1, .1 + legend_space, .8, .8 - legend_space], **axis_settings)
    else:
        ax = fig.get_axes()[-1]
        
    pylabtitle(title_str % locals())

    options = line_settings.copy()
    options.setdefault('linestyle', '-')
    options.setdefault('linewidth', 3)
    options.setdefault('marker', 'None')
    for rxn in reactions:
        options['color'] = colors.next()
        data = rxn[species] * factor
        rxn_sum = rxn.sum()
        reaction_label = rxn.display(digits = None)
        options['label'] = reaction_label
        if rcParams['text.usetex']: options['label'] = '\ce{' + options['label'] + '}'
        ax.plot_date(date_objs, data[slice].repeat(2,0), **options)

    if chem is not None:
        options['color'] = 'black'
        options['marker'] = 'o'
        options['label'] = 'Total Chem'
        try:
            data = array(mech('%s'  % (chem,))[species][slice])
        except:
            warn('Using sum of reactions for %(species)s' % locals())
            reactants = products = []
            if 'p' in species.role() or 'u' in species.role():
                products = species
            if 'u' in species.role() or 'u' in species.role():
                reactants = species
            data = mech.make_net_rxn(reactants, products, logical_and = False, reaction_type = 'kjun')[species][slice]
    
        ax.plot_date(date_objs, data.repeat(2,0) * factor, **options)

    ax.set_xlabel('Time')
    ax.set_ylabel(units)
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
             title = None, filter = True, fig = None, ax = None, secondaryconc = True, cmap = None,
             figure_settings = {}, axis_settings = {},
             line_settings = {'linestyle': '-', 'linewidth': 3, 'marker': 'None'},
             legend_settings = {},
             **kwds):

    """
    mech - permm.Mechanism.Mechanism object
    species - permm.SpeciesGroup.Species object
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

    nlines = len(processes)
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    
    date_objs = get_date_steps(mech.mrg, end_date)
    
    if ax is None:
        if fig is None:
            fig = figure(figsize = (8,6))
            ax = fig.add_axes([0.1, 0.12, 0.8, 0.8], **axis_settings)
        else:
            try:
                ax = fig.get_axes()[0]
            except:
                ax = fig.add_axes([0.1, 0.12, 0.8, 0.8], **axis_settings)
            
    
    if isinstance(secondaryconc, bool):
        if secondaryconc:
            if len(fig.get_axes()) == 2:
                tax = fig.get_axes()[-1]
            else:
                tax = twinx(ax)
        else:
            tax = ax
    else:
        # Assuming secondaryconc is an axis
        tax = secondaryconc
        
    ax.set_title(title_str % locals())
    options = kwds.get('CONC',line_settings.copy())
    if options is not None:
        options.setdefault('color', 'k')
        options.setdefault('label', 'Conc')
        options.setdefault('linestyle', '-')
        options.setdefault('linewidth', 3)
        options.setdefault('marker', 'None')
        init_vals = mech('%s'  % (init,))[species]
        fconc_vals = mech('%s'  % (final,))[species]
        data = concatenate([init_vals[..., :1], fconc_vals[...]], axis = -1) * factor
        tax.plot_date(get_dates(mech.mrg), data, **options)
    units = kwds.get('units', None)
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
        var = mech('(%s)' % (process,))
        units = units or var.get_units(species)
        data = var[species].repeat(2,0) * factor
        if data.nonzero()[0].any() or not filter:
            ax.plot_date(date_objs, data, **options)
            
    if ax.get_xlabel() == '':
        ax.set_xlabel('Time')
    if ax.get_ylabel() == '':
        ax.set_ylabel('processes (%s)' % units)

    if not tax is ax:
        tax.set_ylabel('concentration')

    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    ax.set_xlim(date2num(date_objs[0]), date2num(date_objs[-1]))

    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)

    fig.autofmt_xdate()
    legend_settings.setdefault('prop', FontProperties(size=10))
    
    if tax == ax:
        fig.legend(**legend_settings)
    else:
        legend_lines = [l for l in tax.lines + ax.lines if l.get_label()[:5] != '_line' and l.get_label() is not None]
        legend_labels = [l.get_label() for l in legend_lines]
        fig.legend(legend_lines, legend_labels, **legend_settings)
    return fig

def _add_mech(conf):
    if conf.has_key('mech'):
        mech = conf['mech']
    else:
        from permm import get_pure_mech
        mech = conf['mech'] = get_pure_mech('_'.join([conf['mechanism'].lower(), conf['model'].lower()]))
        if isinstance(conf['mrgfile'], (PseudoNetCDFFile, InstanceType)):
            mrg_file = conf['mrgfile']
        else:
            mrg_file = NetCDFFile(conf['mrgfile'],'r')

        mech.set_mrg(mrg_file)

def phy_plots(conf, filter = True, fmt = 'pdf', fig_append = 'IPR'):
    _add_mech(conf)
    
    mech = conf['mech']
    date_objs = get_date_steps(mech.mrg, conf)

    for species, species_options in conf['species'].iteritems():
        try:
            fig = phy_plot(conf, mech,  date_objs, species, species_options, filter = filter)
            fig.savefig(os.path.join(conf['outdir'], '%s_%s.%s' % (species, fig_append, fmt)), format = fmt)
        except KeyError, detail:
            warn(detail)
    
import unittest

class TestPlots(unittest.TestCase):
    def setUp(self):
        from permm.mechanisms import small_strato
        from PseudoNetCDF import PseudoNetCDFFile
        from numpy import arange, newaxis
        from numpy.random import normal, poisson, random, seed
        
        self.mech = small_strato()
        mrg = self.mrg = PseudoNetCDFFile()
        
        mrg.createDimension('TSTEP', 10)
        mrg.createDimension('ROW', 3)
        mrg.createDimension('LAY', 4)
        mrg.createDimension('COL', 5)
        mrg.createDimension('VAR', 10)
        mrg.createDimension('DATE-TIME', 2)
        
        tflag = mrg.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
        tflag.units = '<YYYYJJJ, HHMMSS>'
        tflag.long_name = tflag.var_desc = 'time '
        
        tflag[:, :, 0] = 2004001
        tflag[:, :, 1] = arange(10)[:, newaxis] * 10000
        
        for i in range(1, 11):
            var = mrg.createVariable('IRR_%d' % i, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'))
            var.units = 'ppt'
            var.long_name = var.var_desc = 'Integrated rate for IRR ordinal %d' % i
            seed(1)
            if i % 2 == 0:
                data = arange(3 * i, 3 * i + 10)[:, newaxis, newaxis, newaxis].repeat(4, 1).repeat(3, 2).repeat(5, 3)
            elif i % 1 == 0:
                data = arange(2 * i, 2 * i + 10)[:, newaxis, newaxis, newaxis].repeat(4, 1).repeat(3, 2).repeat(5, 3)
            else:
                data = arange(1 * i, 1 * i + 10)[:, newaxis, newaxis, newaxis].repeat(4, 1).repeat(3, 2).repeat(5, 3)

            var[:] = data * {True: -1, False: 1}[i % 2 == True]

        self.mech.set_mrg(mrg)

    def testIRR(self):
        from tempfile import TemporaryFile
        tf = TemporaryFile()
        fig = self.mech("plot_rxns(O1D, O1D, False, slice = (slice(None), 0, 1, 2), show = False)")
        fig.savefig(tf)
        tf.seek(0,0)
        tfdata = tf.read()
        testdata = file('testdata/test1.png').read()
        import pdb; pdb.set_trace()
        self.assertTrue(tfdata == testdata)

if __name__ == '__main__':
    unittest.main()