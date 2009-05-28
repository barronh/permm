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
from SpeciesGroup import Species
from ProcessGroup import Process
from ReactionGroup import ReactionFromString
from IPRArray import IPR

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
        for spc in yaml_file.species_list:
            self.species_dict[spc] = Species(name = spc, names = [spc], stoic = [1])
        
        self.reaction_dict = AttrDict()
        for rxn_name, rxn_str in yaml_file.reaction_list.iteritems():
            self.reaction_dict[rxn_name] = ReactionFromString(rxn_str)

        self.__reaction_data = self.reaction_dict

        for grp_name, spc_grp in yaml_file.get('species_group_list',[]):
            self.species_dict[grp_name] = eval(spc_grp,None, self.species_dict)
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
        
    def find_rxns(self, reactants, products, logical_and = True):
        """
        Get reactions where filter is true
        
        When logical_and is true, filter requires reaction to have
        all reactants AND all products.
        
        When logical_and is False, filter requires reaction to have
        all reactants OR all products.
        """
        if isinstance(reactants, (Species, str)):
            reactants = [reactants]
        
        reactants = [self.__ensure_species(spc) for spc in reactants]
            
        if isinstance(products, (Species, str)):
            products = [products]

        products = [self.__ensure_species(spc) for spc in products]
            
        result = [(rxn_name, rxn) for rxn_name, rxn in self.reaction_dict.iteritems()]


        if reactants != []:
            reactant_result = [set([rxn_name for rxn_name, rxn in result if rxn.has_rct(rct)]) for rct in reactants]
            reactant_result = list(reduce(lambda x,y: set.intersection(x,y), reactant_result))
        else:
            reactant_result = [rn for rn, rv in result]
        
        if products != []:
            product_result = [set([rxn_name for rxn_name, rxn in result if rxn.has_prd(prd)]) for prd in products]
            product_result = list(reduce(lambda x,y: set.intersection(x,y), product_result))
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
        
    def make_net_rxn(self, reactants, products, logical_and = True):
        """
        Create and return a net reactions that meet the reactants and
        products filter (see find_rxns)
        """
        rxns = self.find_rxns(reactants, products, logical_and)
        result = self(' + '.join(rxns))
        
        return result

    def print_rxns(self, reactants, products, logical_and = True):
        """
        Print reactions that meet the reactants and
        products filter (see find_rxns)
        """
        for rxn in self.find_rxns(reactants, products, logical_and):
            print rxn, self.reaction_dict[rxn]
        
    def print_nrxns(self, reactants, products, logical_and = True, factor = 1.):
        """
        Print net reactions that meet the reactants and
        products filter (see find_rxns)
        """
        for rxn in self.find_rxns(reactants, products, logical_and):
            print rxn, self.nreaction_dict[rxn].sum() * factor
        
    def set_mrg(self,mrg, use_net_rxns = True):
        """
        Add process analysis from a 1D merged IRR/IPR file
        """
        self.mrg = mrg
        self.set_irr(mrg.variables['IRR'], mrg.Reactions.split())
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
            self.process_dict[prc_name] = eval(prc,{},self.process_dict)
            self.process_dict[prc_name].name = prc_name
            
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
