from utils import AttrDict
from collections import defaultdict
import operator
import yaml
import re
import sys
from numpy import dtype, \
                  array, \
                  ndarray
from warnings import warn
from SpeciesGroup import Species, species_sum
from ProcessGroup import Process
from ReactionGroup import ReactionFromString
from IPRArray import IPR
from graphing.timeseries import irr_plot, phy_plot

__all__ = ['Mechanism']

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
            for spc in yaml_file.species_list:
                self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1])
        
        self.reaction_dict = AttrDict()
        reaction_species = []
        for rxn_name, rxn_str in yaml_file.reaction_list.iteritems():
            rxn = self.reaction_dict[rxn_name] = ReactionFromString(rxn_str)
            reaction_species += rxn.species()
        
        reaction_species = set(reaction_species)
        for spc in [spc for spc in reaction_species if not self.species_dict.has_key(spc)]:
            self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1])

        self.__reaction_data = self.reaction_dict

        for spc_grp_def in yaml_file.get('species_group_list',[]):
            grp_name = spc_grp_def.split('=')[0].strip()
            if (spc_grp_def.count('+') + spc_grp_def.count('-')) > 4:
                spc_grp_def = '=species_sum(['.join(spc_grp_def.replace('+',',').replace('-',',').split('='))+'])'
            exec(spc_grp_def, None, self.species_dict)
            self.species_dict[grp_name].name = grp_name

        self.net_reaction_dict = yaml_file.get('net_reaction_list',{})
            
    def __call__(self,expr):
        """
        Evaluate string expression in the context of mechanism
        species, species groups, reactions, net reactions and processes
        """
        return eval(expr, None, self)
    
    def __getitem__(self,item):
        """
        Provide a single getitem interface for mechanism species, species 
        groups, reactions, net reactions and processes
        """
        try:
            result = eval(item, self.__reaction_data, self.species_dict)
            self.__reaction_data.pop('__builtins__')
        except NameError:
            result = eval(item, self.species_dict, self.__pa_data)
            self.species_dict.pop('__builtins__')
        return result
        
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
        
    def find_rxns(self, reactants = [], products =[], logical_and = True):
        """
        Get reactions that meet the criteria specified by reactants, products and logical_and.
        
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
            result = list(set(reactant_result).intersection(product_result))
        else:
            result = list(set(reactant_result+product_result))

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
        
    def make_net_rxn(self, reactants = [], products = [], logical_and = True):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        rxns = self.find_rxns(reactants, products, logical_and)
        result = self(' + '.join(rxns))
        
        return result

    def print_net_rxn(self, reactants = [], products = [], logical_and = True):
        """
        Sum each reaction in find_rxns(reactants, procucts, logical_and)
        and return a net reaction
        
        for more information see find_rxns
        """
        net_rxn = self.make_net_rxn(reactants, products, logical_and)
        
        print net_rxn

    def print_rxns(self, reactants = [], products = [], logical_and = True):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction
        
        for more information see find_rxns
        """
        for rxn in self.find_rxns(reactants, products, logical_and):
            print rxn, self.reaction_dict[rxn]

    def plot_rxns(self, reactants = [], products = [], logical_and = True, plot_spc = None, path = None, **conf):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        plot the production/consumption rate for plot_spc.  If plot_spc is
        not specified, use the first product provided.  If no products, use 
        the first reactant.
        
        for more information see find_rxns
        """
        def ensure_list(x):
            if isinstance(x, Species):
                return [x]
            else:
                return x
        if plot_spc is None:
            plot_spc = (ensure_list(products)+ensure_list(reactants))[0]
        reactions = self.find_rxns(reactants, products, logical_and)
        
        if reactions == []:
            raise ValueError, "Your query didn't match any reactions; check your query and try again (try print_rxns)."
        fig = irr_plot(self, reactions, plot_spc, **conf)
        if path is not None:
            fig.savefig(path)
        else:
            from pylab import show
            show()
        
        return fig
        
    def print_nrxns(self, reactants = [], products = [], logical_and = True, factor = 1.):
        """
        For each reaction in find_rxns(reactants, procucts, logical_and),
        print the reaction summed for the entire timeseries.
        
        for more information see find_rxns
        """
        if not hasattr(self,'nreaction_dict'):
            raise ValueError, "Net reactions are only available when IRR has been loaded"
        for rxn in self.find_rxns(reactants, products, logical_and):
            print rxn, self.nreaction_dict[rxn].sum() * factor
        
    def set_mrg(self,mrg, use_net_rxns = True, use_irr = True, use_ipr = True):
        """
        Add process analysis from a 1D merged IRR/IPR file
        """
        self.mrg = mrg
        if use_irr:
            self.set_irr(mrg.variables['IRR'], mrg.Reactions.split(), use_net_rxns = use_net_rxns)
        if use_ipr:
            self.set_ipr(mrg.variables['IPR'], mrg.Species.split(), mrg.Process.split())
        
    def set_irr(self,irr, ReactionNames, use_net_rxns = True):
        """
        Add process analysis from a 2D merged IRR array dim(TIME,RXN)
        """
        irr_type = dtype(dict(names = ReactionNames, formats = 'f'*len(ReactionNames)))
        class irr_array(ndarray):
            pass
            
        self.irr = array(irr).view(dtype = irr_type).squeeze().view(type = irr_array)
        self.irr.units = irr.units
        self.nreaction_dict = AttrDict()
        for rxn_name, rxn in self.reaction_dict.iteritems():
            self.nreaction_dict[rxn_name] = rxn * self.irr[rxn_name]

        self.__reaction_data = self.nreaction_dict
        
        if use_net_rxns:
            for nrxn_name, nrxn in self.net_reaction_dict.iteritems():
                self.nreaction_dict[nrxn_name] = eval(nrxn, None, self.nreaction_dict)
        
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
            
        self.ipr = IPR(array(ipr), species, processes)
        self.ipr.units = ipr.units
        
        def pa_dict(item):
            if self.process_dict.has_key(item):
                return self.ipr[self.process_dict[item]]
            else:
                return self.irr[item]
        
        class key_specific_default(defaultdict):
            def __missing__(self, key):
                if self.default_factory is None:
                    raise KeyError(key)
                self[key] = value = self.default_factory(key)
                return value

        self.__pa_data = key_specific_default(pa_dict)

    def add_rxn(self, rxn_key, rxn_str):
        self.reaction_dict[rxn_key] = ReactionFromString(rxn_str)
        if hasattr(self, 'nreaction_dict'):
            self.nreaction_dict[rxn_key] = ReactionFromString(rxn_str)
            self.nreaction_dict[rxn_key] *= self.irr[rxn_key]
    
    def plot_proc(self, species, path = None, **kwds):
        """
        species - perm.SpeciesGroup.Species object
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