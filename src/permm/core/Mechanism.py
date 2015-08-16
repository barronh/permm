import operator
import yaml
import re
import sys
from numpy import * # Explicitly using dtype, array and ndarray; providing all default numpy to __call__interface
from warnings import warn

from PseudoNetCDF.sci_var import PseudoNetCDFVariable

from Species import Species, species_sum
from Reaction import Reaction
from IPRArray import Processes_ProcDelimSpcDict, Process

from permm.graphing.timeseries import irr_plot, phy_plot, plot as tplot
from permm.Shell import load_environ
from permm.netcdf import NetCDFVariable

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
        
        self.__yaml_file = yaml_file
        self.mechanism_comment = yaml_file.get('comment', '')
        self.species_dict = dict()

        for spc, spc_def in yaml_file.get('species_list', {}).iteritems():
            self.species_dict[spc] = Species(spc + ': ' + spc_def)
                        
        self.reaction_dict = dict()
        reaction_species = []
        for rxn_name, rxn_str in yaml_file.get('reaction_list', {}).iteritems():
            rxn = self.reaction_dict[rxn_name] = Reaction(rxn_str)
            reaction_species += rxn.species()
        
        reaction_species = set(reaction_species)
        for spc in [spc for spc in reaction_species if not self.species_dict.has_key(spc)]:
            self.species_dict[spc] = Species(spc + ': IGNORE')

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
        try:
            return self.variables[expr]
        except:
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
        
        rxns = self.find_rxns(reactants = spc)
        return self.__add_spc_to_reactions(rxns, spc)

    def add_prd_to_reactions(self, spc):
        """
        Add spc to reactions where components are products
        see __ensure_species
        """
        spc = self.__ensure_species(spc)
        
        rxns = self.find_rxns(reactants = [],products = spc)
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
    
    def get_rxns(self, reactants = [], products = [], logical_and = True, reaction_type = None):
        """
        Get reaction objects that meet the criteria specified by reactants, products and logical_and.
        
        reactants - filter mechanism reactions for reactions with all reactant(s)
        products - filter mechanism reactions for reactions with all product(s)
        logical_and - a boolean indicating how to combine reactant and product filters
                      True: reaction is in both filters (i.e. reactants AND products)
                      False: reaction is in either filter (i.e. reactants OR products)
        """
        return [self.reaction_dict[rk] for rk in self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)]

    def get_irrs(self, reactants = [], products = [], logical_and = True, reaction_type = None, sortby = None):
        """
        Get reaction objects that meet the criteria specified by reactants, products and logical_and.
        
        reactants - filter mechanism reactions for reactions with all reactant(s)
        products - filter mechanism reactions for reactions with all product(s)
        logical_and - a boolean indicating how to combine reactant and product filters
                      True: reaction is in both filters (i.e. reactants AND products)
                      False: reaction is in either filter (i.e. reactants OR products)
        """
        result = [self(rk) for rk in self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)]
        if sortby is None:
            return result
        else:
            result = [(r[sortby].sum(), r) for r in result]
            result.sort()
            result = [r for s, r in result]
            return result

    def filter_rxns(self, spcs = None, reaction_type = None):
        """
        Return a list of filtered reactions.
        
        
        spc - Species with roles specified; single or list of species with roles
              roles and exclude property affect the outcome.
        reaction_type - some combination of kjn: thermal (k); photolysis (j); net (n)
        
        
        """
        result = [(rxn_name, rxn) for rxn_name, rxn in self.reaction_dict.iteritems()]
        if not reaction_type is None:
            result = [(rn, rx) for rn, rx in result if rx.reaction_type in reaction_type]
        
        if not spcs is None:
            if not isinstance(spcs, list):
                spcs = [spcs]
            spcs = [self.__ensure_species(spc) for spc in spcs]
        
            result = [(rn, rx) for rn, rx in result if all([spc in rx for spc in spcs])]
        return [rn for rn, rx in result]
    
    def find_rxns(self, reactants = [], products = [], logical_and = True, reaction_type = None):
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
            reaction_list = set(reactant_result + product_result)

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
    
    def subst_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None, name = None, netspc = False):
        """
        Substitute net reaction for reaction in mechanism.
        
        1. Find all matching reactions.
        2. Create a net reactions from all matching reactions.
        3. Remove all individual reactions (from e.g., irr_dict, reaction_dict, and variables)
        4. Add net reaction from step 2 (to e.g., irr_dict, reaction_dict, and variables)
        
        Returns: None
        """
        rxns = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
        if len(rxns) > 1:
            eval_str = ' + '.join(rxns)
            if name is None:
                name = eval_str.replace(' + ', 'pl')
            nrxn = eval(eval_str, globals(), self.irr_dict)
            if netspc:
                nrxn = nrxn.net()
            
            for rxn in rxns:
                del self.irr_dict[rxn]
                del self.reaction_dict[rxn]
                del self.variables[rxn]
            self.irr_dict[name] = nrxn
            self.reaction_dict[name] = nrxn.sum()
            load_environ(self, self.variables)

    def subst_these_rxns(self, rxns, name = None):
        """
        
        Synopsis
        Substitute net reaction for reactions in mechanism.
        
        1. Create a net reactions from all reactions (rxns).
        2. Remove all individual reactions (from e.g., irr_dict, reaction_dict, and variables)
        3. Add net reaction from step 1 (to e.g., irr_dict, reaction_dict, and variables)
        
        Requires:
            rxns - iterable of reaction names
            name - key for new reaction
        Returns:
            None
        """
        if all([isinstance(rxn, Reaction) for rxn in rxns]):
            irrs = rxns
            rxns = []
            for irr in irrs:
                for rxnkey, crxn in self.irr_dict.iteritems():
                    if irr is crxn:
                        rxns.append(rxnkey)
                        break
                else:
                    raise KeyError('Could not match %s with a active reaction')
        elif all([isinstance(rxn, str) for rxn in rxns]):
            pass
        else:
            raise ValueError('Reactions must be all strings (i.e. keys) or objects')
                
        if len(rxns) > 1:
            eval_str = ' + '.join(rxns)
            if name is None:
                name = eval_str.replace(' + ', 'pl')
            nrxn = eval(eval_str, globals(), self.irr_dict)
            for rxn in rxns:
                del self.irr_dict[rxn]
                del self.reaction_dict[rxn]
                del self.variables[rxn]
            self.irr_dict[name] = nrxn
            self.reaction_dict[name] = nrxn.sum()
            load_environ(self, self.variables)
        
        
    def make_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        rxns = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
        if len(rxns) == 0:
            return Reaction(stoic = dict())
        result = self(' + '.join(rxns))
        
        return result

    def make_net_filter_rxn(self, spcs = [], reaction_type = None):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        rxns = self.filter_rxns(spcs = spcs, reaction_type = reaction_type)
        if len(rxns) == 0:
            return Reaction(stoic = dict())
        result = self(' + '.join(rxns))
        
        return result

    def print_net_rxn(self, reactants = [], products = [], logical_and = True, reaction_type = None):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        net_rxn = self.make_net_rxn(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
        
        print net_rxn

    def print_rxns(self, reactants = [], products = [], logical_and = True, reaction_type = None, sortby = None, reverse = False, digits = -1, nspc = 10000):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction
        
        for more information see find_rxns
        """
        rxns = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)

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

    def plot_rxn_list(self, reactions, plot_spc = None, path = None, show = False, **kwds):
        """
        Returns:
            a matplotlib figure with a line for each reaction indexed by plot_spc
        
        Parameters:
            reactions - list of reaction objects
            plot_spc - species to plot
            path - path to savefigure
            show - try to display figure interactively

        Other keywords:
            accepts any keyword accepted by permm.graphing.timeseries.irr_plot
        """
        if plot_spc is None:
            plot_spc = self((reactions[0][0].products() + reactions[0][0].reactants())[0])
            
        fig = irr_plot(self, reactions = reactions, species = plot_spc, **kwds)
        if path is not None and show:
            fig.savefig(path)
        else:
            from pylab import show
            show()
        
        return fig

    def plot(self, y, path = None, stepped = True, end_date = True, time_slice = slice(None), figure_settings = {}, axis_settings = {}, line_settings = {}):
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
                  plot_spc = None, combine = [()], nlines = 8, **kwds):
        """
        Creates figure of reactions
        
        Steps:
        1. Query reactions using get_rxns(reactants, procucts, logical_and)
        2. net reaction sets specified in combine (i.e. combine = [('IRR_1', 'IRR_2'), ...])
        3. sort queried reactions by absolute plot_spc change (+/-)
        3. combine all reactions that are not in the top (nlines - 1)
        4. plot the plot_spc change from each (nlines - 1) individual reactions
           and the 1 net reaction
        
        Returns:
            matplotlib figure with each queried reaction plotted for production/consumption
            of plot_spc.
        
        Parameters:
            reactants - instance or list of Species objects to require
            products - instance or list of Species objects to require
            logical_and - if true, require reactants and products; false, require reactants or products
            reaction_type - some combination of kjn: thermal (k); photolysis (j); net (n)
            plot_spc - Species instance to be plotted, defaults to first product or reactant
            combine - reaction keys to combine into net reactions
            nlines - maximum reactions to plot; selected based on maximum absolute plot_spc change
        
        Other keywords:
            permm.Mechanism.plot_rxn_list or 
            permm.graphing.timeseries.irr_plot
        
        """
        if plot_spc is None:
            plot_spc = (_ensure_list(products)+_ensure_list(reactants))[0]
        reactions = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
        
        if reactions == []:
            raise ValueError, "Your query didn't match any reactions; check your query and try again (try print_rxns)."

        reactions1 = reactions
        reactions = [ rxn for rxn in reactions if rxn not in reduce(operator.add, combine) ]
        nlines = min(nlines, len(reactions)+1)
        if combine != [()] and reactions != reactions1:
            reactions = reactions + map(lambda t2: '+'.join(t2), combine)
        aslice = kwds.get('slice', slice(None))
        reactions = [ (abs(self('(%s)' % (rxn))[aslice].get_spc(plot_spc, float64(0.))).sum(),rxn) for rxn in reactions]
    
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
        
        return self.plot_rxn_list(reactions = reactions, plot_spc = plot_spc, **kwds)
        
    def print_irrs(self, reactants = [], products = [], logical_and = True, reaction_type = None, factor = 1., sortby = None, reverse = False, slice = slice(None), nspc = 100000, digits = -1, formatter = 'g'):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction summed for the entire timeseries.
        
        for more information see find_rxns
        """
        if not hasattr(self,'irr_dict'):
            raise ValueError, "Net reactions are only available when IRR has been loaded"
            
        rxns = self.find_rxns(reactants = reactants, products = products, logical_and = logical_and, reaction_type = reaction_type)
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
        
    def set_mrg(self, mrg, use_net_rxns = True, use_irr = True, use_ipr = True):
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
                self.set_ipr(mrg.variables['IPR'])
            except:
                self.set_ipr()
                        
        
        load_environ(self, self.variables)
        
    def set_irr(self, irr = None, ReactionNames = None, use_net_rxns = True):
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
        
    def set_ipr(self, ipr = None, processes = None):
        """
        Add process analysis from a 3D merged IPR array (TIME,SPC,PROC)
        """
        if ipr is None:
            if hasattr(self.mrg, 'Processes'):
                processes = self.mrg.Processes.split()
                self.process_dict = Processes_ProcDelimSpcDict(processes, self.mrg.variables)
            else:
                warn("Unable to load IPR; missing Species and Processes attributes")
                return
        elif isinstance(ipr, dict):
            if processes is None:
                raise ValueError, "When ipr is a dictionary, processes must be provided as a list of process names"
            self.process_dict = Processes_ProcDelimSpcDict(processes, self.mrg.variables)
        elif isinstance(ipr, (PseudoNetCDFVariable, NetCDFVariable)):
            self.process_dict = {}
            for pi, prc in enumerate(processes):
                self.process_dict[prc] = Process(prc, default_unit = getattr(ipr, 'units', 'Unknown'), **dict([(spc, ipr[si, pi]) for si, spc in enumerate(species)]))
        else:
            return

        for prc_name, prc in self.__yaml_file.get('process_group_list', {}).iteritems():
            try:
                self.process_dict[prc_name] = eval(prc,{},self.process_dict)
                self.process_dict[prc_name].name = prc_name
            except:
                warn("Cannot create %s process group" % prc_name)
            
        # Add extra species for names in IPR
        spcs = set(reduce(list.__add__, [[spc for spc in proc.keys()] for proc in self.process_dict.values()]))
        new_spcs = spcs.difference(self.species_dict.keys())
        for name in new_spcs:
            if not self.species_dict.has_key(name):
                self.species_dict[name] = Species(name + ': IGNORE')

    def add_rxn(self, rxn_key, rxn_str):
        """
        Synopsis:
            add reaction with name rxn_key and definition rxns_str
        
        Requires:
            rxn_key - string to be used as key for reaction.
            rxn_str - string that defines the reaction (e.g., N2O + 2 H2O ->[k] HNO3
        """
        self.reaction_dict[rxn_key] = Reaction(rxn_str)
        if hasattr(self, 'irr_dict'):
            self.irr_dict[rxn_key] = Reaction(rxn_str)
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
            if kwds.get('show', True):
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
        