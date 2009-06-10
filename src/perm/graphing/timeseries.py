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
from numpy import array, concatenate, zeros, arange
from pylab import figure, title, plot_date, savefig, legend, axis, grid, axes, xlabel, ylabel, subplot
from matplotlib.dates import DateFormatter, date2num
from matplotlib.cm import get_cmap
from matplotlib.font_manager import FontProperties
import re
import operator
import os

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
    date_objs = get_date(mech.mrg, conf)
    units = mech.irr.units
    nlines = conf.get('nlines',8)
    combine = conf.get('combine',[()])
    fig = conf.get('fig', None)
    cmap = conf.get('cmap', None)
    title_str = conf.get('title', 'Plot of %s for %d Reactions' % (species.name, len(reactions)))
    chem = conf.get('chem', 'Chemistry')
    
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    if fig is None:
        fig = figure()
    ax = axes([0.1,0.4,.8,.5], **conf)
    grid(True)
    title(title_str % locals())

    reactions = [ rxn for rxn in reactions if rxn not in reduce(operator.add, combine) ]
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
        data = mech('(%s)' % (rxn, ))[species]
        options['label'] = re.compile('\d+.\d+\*').sub('', str(mech(rxn).sum())).replace(' ', '')
        plot_date(date_objs, data.repeat(2,0), **options)

    data = zeros(data.shape, dtype = data.dtype)
    for rxn in reactions[nlines-1:]:
        data += mech('(%s)' % (rxn,))[species]
    options['label'] = 'other'
    plot_date(date_objs, data.repeat(2,0), **options)
    
    options['color'] = 'black'
    options['marker'] = 'o'
    options['label'] = 'Chem'
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
    
    legend(loc=(0,-0.8), prop = FontProperties(size=10))
    return fig

def phy_plot(mech, species, **conf):
    """
    conf - configuration obect that has title, 
           init, final, process, and species
           * title - string template that has access
                     to all local variables by name
           * mech - net_balance.Mechanism.Mechanism object
           * species - sp
    """
    filter = conf.get('filter', True)
    fig = conf.get('fig', None)
    init = conf.get('init', 'Initial')
    final = conf.get('final', 'Final')
    exclude = conf.get('exclude', 'UCNV AVOL DATE TIME I J K'.split())
    title_str = conf.get('title', '%s Process' % species.name)
    processes = conf.get('process', dict([(p,{}) for p in set(mech.ipr.dtype.fields.values()[0][0].fields.keys()).difference([init, final]+exclude)]))
    cmap = conf.get('cmap', None)
    units = mech.ipr.units
    nlines = len(processes)
    colors = iter(get_cmap(cmap)(arange(nlines, dtype = 'f')/(nlines-1)))
    
    date_objs = get_date(mech.mrg, conf)
    if fig is None:
        fig = figure()
    
    ax = axes(**conf)

    grid(True)
    title(title_str % locals())
    options = {'color': 'k'}
    options.setdefault('linestyle', '-')
    options.setdefault('linewidth', 3)
    options.setdefault('marker', 'o')
    data = mech('%s'  % (init,))[species].array()
    plot_date(date_objs[::2], data, **options)
    options.setdefault('label', 'Conc')
    options['marker'] = 'x'
    data = mech('%s'  % (final,))[species].array()
    plot_date(date_objs[1::2], data, **options)
    
    for process in processes.keys():
        options = processes.get(process, {})
        options.setdefault('color', colors.next())
        options.setdefault('linestyle', '-')
        options.setdefault('linewidth', 3)
        options.setdefault('label', process)
        options.setdefault('marker', 'None')
        data = mech('(%s)' % (process,))[species].array().repeat(2,0)
        if data.nonzero()[0].any() or not filter:
            plot_date(date_objs, data, **options)
            
    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%jT%H'))
    ax.set_xlim(date2num(date_objs[0]), date2num(date_objs[-1]))
    if conf.has_key('ylim'):
        ax.set_ylim(*conf['ylim'])

    fig.autofmt_xdate()
    legend()
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
