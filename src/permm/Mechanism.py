from utils import AttrDict
from collections import defaultdict
import operator
import yaml
import re
import sys
from numpy import * # Explicitly using dtype, array and ndarray; providing all default numpy to __call__interface
from warnings import warn
from SpeciesGroup import Species, species_sum
from ProcessGroup import Process
from ReactionGroup import ReactionFromString
from IPRArray import IPR
from graphing.timeseries import irr_plot, phy_plot, plot as tplot
from Shell import load_environ
from PseudoNetCDF.sci_var import PseudoNetCDFVariable
__all__ = ['Mechanism']

_spc_def_re = re.compile(r'(?P<stoic>[-+]?[0-9]*\.?[0-9]+)(?P<atom>\S+)(?=\s*\+\s*)?')
_numre = re.compile('(\d+)')

class Mechanism(object):
    """
    The mechanism object is the central repository for information
    about a chemcial model.  It is aware of species, reactions, species
    groups (e.g. NOy), net reactions (e.g. PAN <-> C2O3 + NO2).  
    
    The mechanism object also has introspective tools.  In addition, it
    provides an interpreter through its call interface.
    
    Reaction data can be augmented with IRR data via set_irr or set_mrg
    """
    def __init__(self, yaml_path):
        """
        Initialization requires a yaml_path or a dictionary.
        The dictionary should have the following keys:
        
        species_list - a list of string names for chemical species
                        (i.e. O3)
        reaction_list - a list of strings representing reactions 
                        (see ReactionGroup.ReactionFromString)
        net_reaction_list - a dictionary of reaction additions 
                            (i.e. NAME: RXN_01 + RXN_02)
        process_group_list - a dictionary of process additions
                            (i.e. VertAdv = Top_Adv + Bottom_Adv)
        """
        import os
        if isinstance(yaml_path,str):
            if os.path.exists(yaml_path):
                yaml_file = yaml.load(file(yaml_path))
            else:
                yaml_file = yaml.load(yaml_path)
        elif isinstance(yaml_path,dict):
            yaml_file = yaml_path
        
        yaml_file = AttrDict(yaml_file)
        self.__yaml_file = yaml_file
        self.species_dict = AttrDict()
        if yaml_file.has_key('species_list'):
            for spc, spc_def in yaml_file.species_list.iteritems():
                if spc_def == 'IGNORE':
                    self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1])
                else:
                    spc_atoms = dict([(atom, eval(stoic)) for stoic, atom in _spc_def_re.findall(spc_def)])
                    self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1], atom_dict = {spc: spc_atoms})
                        
        self.reaction_dict = AttrDict()
        reaction_species = []
        for rxn_name, rxn_str in yaml_file.reaction_list.iteritems():
            rxn = self.reaction_dict[rxn_name] = ReactionFromString(rxn_str)
            reaction_species += rxn.species()
        
        reaction_species = set(reaction_species)
        for spc in [spc for spc in reaction_species if not self.species_dict.has_key(spc)]:
            self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1])

        for spc_grp_def in yaml_file.get('species_group_list',[]):
            grp_name = spc_grp_def.split('=')[0].strip()
            if (spc_grp_def.count('+') + spc_grp_def.count('-')) > 4:
                spc_grp_def = '=species_sum(['.join(spc_grp_def.replace('+',',').replace('-',',').split('='))+'])'
            exec(spc_grp_def, None, self.species_dict)
            self.species_dict[grp_name].name = grp_name

        self.net_reaction_dict = yaml_file.get('net_reaction_list',{})
        self.variables = {}
        load_environ(self, self.variables)
            
    def __call__(self, expr, env = globals()):
        """
        Evaluate string expression in the context of mechanism
        species, species groups, reactions, net reactions and processes
        """
        return eval(expr, env, self.variables)
    
    def __getitem__(self,item):
        """
        Provide a single getitem interface for mechanism species, species 
        groups, reactions, net reactions and processes
        """
        return self.variables[item]
        
    def __add_spc_to_reactions(self,rxn_list, spc):
        """
        Add a species to a reaction as an accumulation of components
          if RXN1 = CH3C(O)CH3 -j> CH3C(O)O2 + CH3O2
          and SPC = Species(CH3C(O)O2 + CH3O2, name = 'AO2')
          
          mech.__add_spc_to_reactions([RXN1],SPC)
          
          would make RXN1 = CH3C(O)CH3 -j> CH3C(O)O2 + CH3O2 + 2*AO2
        
        For more detail, see ReactionGroup.Reaction.__add__ interface
        """
        if not self.species_dict.has_key(spc.name):
            self.species_dict[spc.name] = spc
            
        for rxn in rxn_list:
            self.reaction_dict[rxn] = self.reaction_dict[rxn] + spc
        
        return len(rxn_list)
    
    def __ensure_species(self,spc):
        """
        Return species (created if necessary)
        """
        if isinstance(spc,str):
            spc = self.species_dict[spc]
        return spc
        
    def add_rct_to_reactions(self, spc):
        """
        Add spc to reactions where components are reactants
        see __ensure_species
        """
        spc = self.__ensure_species(spc)
        
        rxns = self.find_rxns(spc,[])
        return self.__add_spc_to_reactions(rxns, spc)

    def add_prd_to_reactions(self, spc):
        """
        Add spc to reactions where components are products
        see __ensure_species
        """
        spc = self.__ensure_species(spc)
        
        rxns = self.find_rxns([],spc)
        return self.__add_spc_to_reactions(rxns, spc)

    def add_spc_to_reactions(self, spc):
        """
        Add spc to reactions where components are reactants
        or products see __ensure_species
        """
        nrcts = self.add_rct_to_reactions(spc)
        nprds = self.add_prd_to_reactions(spc)
        return nrcts+nprds

    def add_spc_to_reactions_with_condition(self, spc, condition):
        """
        spc - must be a species object or string
        condition - must be a function that, given a reaction, evaluates to true or false
        """
        spc = self.__ensure_species(spc)
        
        new_rxn_def = {}
        for rxn_name, rxn in self.reaction_dict.iteritems():
            if condition(rxn):
                rxn = rxn + spc
                new_rxn_def[rxn_name] = rxn
        
        self.reaction_dict.update(new_rxn_def)

        return len(new_rxn_def)
    
    def get_rxns(self, reactants = [], products =[], logical_and = True, reaction_type = None):
        """
        Get reaction objects that meet the criteria specified by reactants, products and logical_and.
        
        reactants - filter mechanism reactions for reactions with all reactant(s)
        products - filter mechanism reactions for reactions with all product(s)
        logical_and - a boolean indicating how to combine reactant and product filters
                      True: reaction is in both filters (i.e. reactants AND products)
                      False: reaction is in either filter (i.e. reactants OR products)
        """
        return [self.reaction_dict[rk] for rk in self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)]

    def get_irrs(self, reactants = [], products =[], logical_and = True, reaction_type = None):
        """
        Get reaction objects that meet the criteria specified by reactants, products and logical_and.
        
        reactants - filter mechanism reactions for reactions with all reactant(s)
        products - filter mechanism reactions for reactions with all product(s)
        logical_and - a boolean indicating how to combine reactant and product filters
                      True: reaction is in both filters (i.e. reactants AND products)
                      False: reaction is in either filter (i.e. reactants OR products)
        """
        return [self(rk) for rk in self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)]

    def find_rxns(self, reactants = [], products =[], logical_and = True, reaction_type = None):
        """
        Get reaction names that meet the criteria specified by reactants, products and logical_and.
        
        reactants - filter mechanism reactions for reactions with all reactant(s)
        products - filter mechanism reactions for reactions with all product(s)
        logical_and - a boolean indicating how to combine reactant and product filters
                      True: reaction is in both filters (i.e. reactants AND products)
                      False: reaction is in either filter (i.e. reactants OR products)
        """
        if isinstance(reactants, (Species, str)):
            reactants = [reactants]
        
        reactants = [self.__ensure_species(spc) for spc in reactants]
            
        if isinstance(products, (Species, str)):
            products = [products]

        products = [self.__ensure_species(spc) for spc in products]
            
        result = [(rxn_name, rxn) for rxn_name, rxn in self.reaction_dict.iteritems()]


        if reactants != []:
            reactant_result = [rxn_name for rxn_name, rxn in result if all([rxn.has_rct(rct) for rct in reactants])]
        else:
            reactant_result = [rn for rn, rv in result]
        
        if products != []:
            product_result = [rxn_name for rxn_name, rxn in result if all([rxn.has_prd(prd) for prd in products])]
        else:
            product_result = [rn for rn, rv in result]
        
            
        if logical_and:
            reaction_list = set(reactant_result).intersection(product_result)
        else:
            reaction_list = set(reactant_result+product_result)

        if reaction_type is not None:
            reaction_result = [rn for rn, rv in result if rv.reaction_type in reaction_type]
            reaction_list = reaction_list.intersection(reaction_result)

        
        result = list(reaction_list)
        result.sort()
        
        return result
    
    def yaml_net_rxn(self, rxns):
        """
        Create the YAML representation of a net reaction for the supplied
        reactions
        """
        species = list(set(reduce(operator.add, [self.reaction_dict[rxn].species() for rxn in rxns])))
        
        for role in ('r','p'):
            for spc in species:
                for rxn in rxns:
                    method = dict( r = 'has_rct', p = 'has_prd')
                    if getattr(self.reaction_dict[rxn], method)(spc):
                        pass
    
    def subst_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None, name = None):
        rxns = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
        if len(rxns) > 1:
            eval_str = ' + '.join(rxns)
            if name is None:
                name = eval_str
            nrxn = eval(eval_str, globals(), self.irr_dict)
            for rxn in rxns:
                del self.irr_dict[rxn]
                del self.reaction_dict[rxn]
            self.irr_dict[name] = nrxn
            self.reaction_dict[name] = nrxn.sum()
            load_environ(self, self.variables)
        
        
    def make_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        rxns = self.find_rxns(reactants, products, logical_and, reaction_type)
        result = self(' + '.join(rxns))
        
        return result

    def print_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        net_rxn = self.make_net_rxn(reactants, products, logical_and, reaction_type)
        
        print net_rxn

    def print_rxns(self, reactants = [], products = [], logical_and = True, reaction_type = None, sortby = None, reverse = False, digits = -1, nspc = 10000):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction
        
        for more information see find_rxns
        """
        rxns = self.find_rxns(reactants, products, logical_and, reaction_type)

        rxns = [(rxn, self.reaction_dict[rxn]) for rxn in rxns]
        try:
            if sortby is None:
                sortby = ([p for p in _ensure_list(products) if not p.exclude]+[r for r in _ensure_list(reactants)  if not r.exclude])[0]
            rxns = [(rxno[sortby], rxn, rxno) for rxn, rxno in rxns]
            rxns.sort(reverse = reverse)
        except:
            rxns = [(_rxn_ordinal(rxn), rxn, rxno) for rxn, rxno in rxns]
            rxns.sort()
            warn("Not all reactions contain %s; check query and/or explicitly define sortby species" % str(sortby))

        for mass, rxn, rxno in rxns:
            print rxn, rxno.display(nspc = nspc, digits = digits)

    def plot_rxn_list(self, reactions, plot_spc = None, path = None,
                      factor = 1, units = None, end_date = True,
                      chem = 'CHEM', title = None, ylim = None, xlim = None,
                      figure_settings = {}, axis_settings = {}, line_settings = {},
                      fig = None, cmap = None, ncol = 2
                    ):
        if plot_spc is None:
            plot_spc = self((reactions[0][0].products() + reactions[0][0].reactants())[0])
            
        fig = irr_plot(self, reactions = reactions, species = plot_spc,
                       factor = factor, units = units, end_date = end_date,
                       chem = chem, title = title, ylim = ylim, xlim = xlim,
                       figure_settings = figure_settings, axis_settings = axis_settings, line_settings = line_settings,
                       fig = fig, cmap = cmap, ncol = ncol
                      )
        if path is not None:
            fig.savefig(path)
        else:
            from pylab import show
            show()
        
        return fig

    def plot(self, y, path = None, stepped = True, end_date = True, time_slice = None, figure_settings = {}, axis_settings = {}, line_settings = {}):
        """
        plot creates a timeseries plot of any numerical array.
        """
        fig = tplot(self, y, stepped = stepped, end_date = end_date,
                    figure_settings = figure_settings,
                    axis_settings = axis_settings,
                    line_settings = line_settings)
        if path is not None:
            fig.savefig(path)

        return fig

    def plot_rxns(self, reactants = [], products = [], logical_and = True, reaction_type = None,
                  plot_spc = None, path = None,
                  combine = [()], nlines = 8,
                  factor = 1, units = None, end_date = True,
                  chem = 'CHEM', title = None, ylim = None, xlim = None, 
                  figure_settings = {}, axis_settings = {}, line_settings = {},
                  fig = None, cmap = None, ncol = 2
                  ):
        """
        For the top nlines reactions in find_rxns(reactants, procucts, logical_and),
        plot the production/consumption rate for plot_spc.  The top lines are selected
        based on total mass throughput of the species being plotted (i.e. plot_spc).
        If plot_spc is not specified, use the first product provided.  If no products, 
        use the first reactant.
        
        For each set of reactions in combine (i.e. combine = [('IRR_1', 'IRR_2'), ...])
        create a net reaction and plot that instead.  For each reaction, over nlines
        create 
        
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
        if plot_spc is None:
            plot_spc = (_ensure_list(products)+_ensure_list(reactants))[0]
        reactions = self.find_rxns(reactants, products, logical_and, reaction_type)
        
        if reactions == []:
            raise ValueError, "Your query didn't match any reactions; check your query and try again (try print_rxns)."

        reactions = [ rxn for rxn in reactions if rxn not in reduce(operator.add, combine) ]
        nlines = min(nlines, len(reactions)+1)
        if combine != [()]:
            reactions = reactions + map(lambda t2: '+'.join(t2), combine)
        
        reactions = [ (abs(self('(%s)' % (rxn))[plot_spc]).sum(),rxn) for rxn in reactions]
    
        reactions.sort(reverse = True)

        reactions = [r for v,r in reactions]
        for rxn in reactions[nlines-1:]:
            try:
                other = other + self('(%s)' % (rxn,))
            except:
                other = self('(%s)' % (rxn,))
        
        try:
            reactions = [self('(%s)' % (rxn, )) for rxn in reactions[:nlines-1]] + [other]
        except:
            reactions = [self('(%s)' % (rxn, )) for rxn in reactions[:nlines-1]]
        
        return self.plot_rxn_list(reactions = reactions, plot_spc = plot_spc, path = path,
                                  factor = factor, units = units, end_date = end_date,
                                  chem = chem, title = title, ylim = ylim, xlim = xlim, 
                                  figure_settings = figure_settings, axis_settings = axis_settings,
                                  line_settings = line_settings, fig = fig, cmap = cmap, ncol = ncol
                                 )
        
    def print_irrs(self, reactants = [], products = [], logical_and = True, reaction_type = None, factor = 1., sortby = None, reverse = False, slice = None, nspc = 100000, digits = -1, formatter = 'g'):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction summed for the entire timeseries.
        
        for more information see find_rxns
        """
        if not hasattr(self,'irr_dict'):
            raise ValueError, "Net reactions are only available when IRR has been loaded"
            
        rxns = self.find_rxns(reactants, products, logical_and, reaction_type)
        irrs = [(rxn, self.irr_dict[rxn][slice].sum()) for rxn in rxns]
        
        try:
            if sortby is None:
                sortby = ([p for p in _ensure_list(products) if not p.exclude]+[r for r in _ensure_list(reactants)  if not r.exclude])[0]
            irrs = [(irr[sortby], rxn, irr) for rxn, irr in irrs]
            irrs.sort(reverse = reverse)
        except:
            irrs = [(_rxn_ordinal(rxn), rxn, irr) for rxn, irr in irrs]
            irrs.sort()
            warn("Not all reactions contain %s; check query and/or explicitly define sortby species" % str(sortby))

        irrs = [(rxn, irr) for mass, rxn, irr in irrs]


        for rxn, irr in irrs:
            print rxn, (irr * factor).display(nspc = nspc, digits = digits)
        
    def set_mrg(self,mrg, use_net_rxns = True, use_irr = True, use_ipr = True):
        """
        Add process analysis from a 1D merged IRR/IPR file
        """
        self.mrg = mrg
        if use_irr:
            try:
                self.set_irr(mrg.variables['IRR'], mrg.Reactions.split(), use_net_rxns = use_net_rxns)
            except:
                self.set_irr()
        if use_ipr:
            try:
                self.set_ipr(mrg.variables['IPR'], mrg.Species.split(), mrg.Process.split())
            except:
                if hasattr(mrg, 'Species'):
                    spcs = mrg.Species.split()
                    prcs = mrg.Processes.split()
                    temp_var = mrg.variables['%s_%s' % (prcs[0], spcs[0])]
                    ipr = PseudoNetCDFVariable(None, 'ipr', temp_var.dtype.char, ('TSTEP', 'SPECIES', 'PROCESS'), units = temp_var.units, values = zeros(temp_var.shape + (len(spcs), len(prcs)), temp_var.dtype))
                    for si, spc in enumerate(spcs):
                        for pi, prc in enumerate(prcs):
                            ipr[..., si, pi] = mrg.variables['%s_%s' % (prc, spc)][:]
    
                    self.set_ipr(ipr, mrg.Species.split(), mrg.Processes.split())
                else:
                    warn("Unable to load IPR; missing Species and Processes attributes")
                        
        
        load_environ(self, self.variables)
        
    def set_irr(self,irr = None, ReactionNames = None, use_net_rxns = True):
        """
        Add process analysis from a 2D merged IRR array dim(TIME,RXN)
        """
        if not irr is None:
            irr_type = dtype(dict(names = ReactionNames, formats = irr[:].dtype.char*len(ReactionNames)))
                
            self.irr = irr[:].view(dtype = irr_type).squeeze().view(type = PseudoNetCDFVariable)
            self.irr.units = irr.units

        self.__use_net_rxns = use_net_rxns
        self.apply_irr()

    def apply_irr(self):
        self.irr_dict = {}
        for rxn_name, rxn in self.reaction_dict.iteritems():
            if hasattr(self, 'irr'):
                try:
                    self.irr_dict[rxn_name] = rxn * self.irr[rxn_name]
                except ValueError, (e):
                    self.irr_dict[rxn_name] = rxn * zeros(self.irr.shape, 'f')
                    warn("IRR does not contain %s: skipped." % rxn_name)
            else:
                try:
                    self.irr_dict[rxn_name] = rxn * self.mrg.variables[rxn_name][:].view(type = PseudoNetCDFVariable)
                except (KeyError, ValueError), (e):
                    warn("IRR does not contain %s: skipped." % rxn_name)
                
        if self.__use_net_rxns and len(self.irr_dict)>0:
            self.nreaction_dict = {}
            for nrxn_name, nrxn in self.net_reaction_dict.iteritems():
                try:
                    self.nreaction_dict[nrxn_name] = eval(nrxn, None, self.irr_dict)
                except Exception, (e):
                    warn("Predefined net rxn %s is not available; %s" % (nrxn_name, str(e)))

        load_environ(self, self.variables)
        
    def set_ipr(self,ipr, species, processes):
        """
        Add process analysis from a 3D merged IPR array (TIME,SPC,PROC)
        """
        self.process_dict={}
        for proc in processes:
            self.process_dict[proc] = Process(name = proc, names = [proc])
        
        for prc_name, prc in self.__yaml_file.get('process_group_list', {}).iteritems():
            try:
                self.process_dict[prc_name] = eval(prc,{},self.process_dict)
                self.process_dict[prc_name].name = prc_name
            except:
                warn("Cannot create %s process group" % prc_name)
            
        self.ipr = IPR(ipr, species, processes, ipr.units)
        self.ipr.units = ipr.units
        
        # Add extra species for names in IPR
        for name in self.ipr.dtype.names:
            if not self.species_dict.has_key(name):
                self.species_dict[name] = Species(name = name, names = [name], stoic = [1.])

    def add_rxn(self, rxn_key, rxn_str):
        self.reaction_dict[rxn_key] = ReactionFromString(rxn_str)
        if hasattr(self, 'irr_dict'):
            self.irr_dict[rxn_key] = ReactionFromString(rxn_str)
            self.irr_dict[rxn_key] *= self.irr[rxn_key]
        load_environ(self, self.variables)
    
    def plot_proc(self, species, path = None, **kwds):
        """
        species - permm.SpeciesGroup.Species object
        path - path for saved figure
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
        fig = phy_plot(self, species, **kwds)
        if path is not None:
            fig.savefig(path)
        else:
            from pylab import show
            show()
        return fig

    def globalize(self, env):
        load_environ(self, env)

def _ensure_list(x):
    if isinstance(x, Species):
        return [x]
    else:
        return x

def _rxn_ordinal(rxnlabel):
    result = _numre.search(rxnlabel)
    if result is None:
        return None
    else:
        val = result.groups()[0]
        try:
            return eval(val)
        except:
            return int(val)
        