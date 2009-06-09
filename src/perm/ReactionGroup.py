from numpy import ndarray, \
                  array, \
                  newaxis, \
                  float64

from utils import AttrDict
from SpeciesGroup import Species
from warnings import warn
import operator
import re

ReactionGroup = str

__all__ = ['Stoic', 'ReactionFromString', 'Reaction', 'ReactionArray']

class Stoic(float64):
    """
    Stoic is a sub-class of float with a role
    property.  The role property is 'r', 'p', or 'u'.
    r = reactant
    p = product
    u = undetermined
    
    Undetermined is useful when a species is part of a 
    net reaction and may play either role
    
    Stoic are created by
        Stoic(2., 'r')
    
    Stoic supports:
        __mul__: Stoic*(Stoic|Number)
        __rmul__: (Stoic|Number)*Stoic
        __add__: Stoic+(Stoic|Number)
    """
    def __new__(subtype,f,role):
        result = float64.__new__(subtype,f)
        result.role = role
        return result

    def __mul__(self,y):
        return Stoic(float64.__mul__(self,y),self.role)

    def __rmul__(self,y):
        return Stoic(float64.__mul__(self,y),self.role)

    def __add__(self,y):
        if isinstance(y,Stoic):
            if y.role == self.role:
                role = self.role
            else:
                role = 'u'
        else:
            role = self.role

        return Stoic(float64.__add__(self,y),role)

def ReactionFromString(rxn_str):
    """
    ReactionFromString is a convenience function.  It creates
    Reaction objects from a string with the following pattern:
        "(?P<reactants>.*)=(?P<rxn_type>[kj])[>]\s*(?P<products>.*)"
    where each reactant or product matches the following pattern:
        "(\s?(?P<sign>[+-])?\s?)?((?P<stoic>\d(\.(\d{1,3}(E\d{2})?)?)?)\*)?(?P<name>[A-Z]\w{0,3})\s?"

    For example:
        OH + OLE =k> 0.8*FORM + 0.33*ALD2 + 0.62*ALDX + 0.8*XO2 + 0.95*HO2 - 0.7 PAR
    """
    result = AttrDict()
    
    reaction_re = re.compile("(?P<reactants>.*)=(?P<rxn_type>[kj])[>]\s*(?P<products>.*)")

    species_re = re.compile("(\s?(?P<sign>[+-])?\s?)?((?P<stoic>\d(\.(\d{1,3}(E\d{2})?)?)?)\*)?(?P<name>[xyA-Z]\w*)(?:[ +=]|$)+",re.M)
    
    reaction_match = reaction_re.match(rxn_str)
    if reaction_match is None:
        raise SyntaxError, "Reactions must match the following patter\n\n<%d> stoic*spc ([+-] stoic*spc)* =[kju]= [stoic*spc] ([+-] stoic*spc)*"

    reactants = reaction_match.groupdict()['reactants']
    reaction_type = reaction_match.groupdict()['rxn_type']
    products = reaction_match.groupdict()['products']

    result.reaction_type = reaction_type
    for spc in species_re.finditer(reactants):
        spc_g = spc.groupdict()
        name = spc_g['name']
        sign = spc_g['sign']
        stoic = spc_g['stoic']
        if sign is None:
            sign = ''

        if stoic is None:
            stoic = '1'
        
        result[name] = Stoic(-float(sign + stoic) + result.get(name,0), 'r')
    
    for spc in species_re.finditer(products):
        spc_g = spc.groupdict()
        name = spc_g['name']
        sign = spc_g['sign']
        stoic = spc_g['stoic']
        if sign is None:
            sign = '+'

        if stoic is None:
            stoic = '1'
        value = float(sign + stoic) + result.get(name, 0)
        role = result.get(name,Stoic(0,'p')).role
        if role == 'r':
            role = 'u'
        result[name] = Stoic(value, role)
        
    return Reaction(result)

class Reaction(AttrDict):
    """
    Reaction is an object that represents reaction groups.  The simplest
    case being a "single reaction" reaction group.
    
    Reaction groups support the following interfaces

        1) indexing (__getitem__) for stoiciometry
        2) multiplication (__mul__ and __rmul__) by arrays and numbers
        3) addition (__add__) of reactions or species
    
    A Reaction can also determine when a species is a reactant, product 
    or unspecified.  When a species is only present as a one role or the other 
    (reactant or product), it is always that role.  When a species is present 
    in the group as both, its current role is determined by the stoichiometry.
    
    Convenience functions:
        has_spc = spc in Reaction.species()
        has_rct = spc in Reaction.reactants()
        has_prd = spc in Reaction.products()
    
    Other functions:
        add_rct_spc
        add_prd_spc
    """
    def __new__(subtype, *args, **kwds):        
        result = AttrDict.__new__(subtype)
        return result
    
    def __init__(self, *args, **kwds):
        AttrDict.__init__(self, *args, **kwds)
        
        has_rxn_type = self.has_key('reaction_type')
        check = reduce(operator.and_,[isinstance(v,Stoic) for k, v in self.iteritems() if k != 'reaction_type'])
        
        if not check or not has_rxn_type:
            raise ValueError, "Reactions require a reaction_type and species keys with stoic"
    
    def __getitem__(self, item):
        if not isinstance(item, Species):
            return AttrDict.__getitem__(self,item)
        elif self.has_key(item.name):
            return AttrDict.__getitem__(self,item.name)
        else:
            species = [name for name in item.names() if name in self.keys()]
            value = sum([item[spc][0] * AttrDict.__getitem__(self, spc) for spc in species])

            first_spc_role = self[species[0]].role
            same_role = array([first_spc_role == self[spc].role for spc in species]).all()
            
            if same_role:
                role = first_spc_role
            else:
                role = 'u'
            return Stoic(value, role = role)
    
    def __str__(self):
        reactants = [(self[rct], rct) for rct in self.reactants()]
        reactants.sort(reverse=True)
        reactants = ' + '.join(['%.5f*%s' % (-1*stoic,rct) for stoic, rct in reactants])

        products = [(self[prd], prd) for prd in self.products()]
        products.sort(reverse=False)
        products = ' + '.join(['%.5f*%s' % (stoic,prd) for stoic, prd in products])


        result = '%s =%s> %s' % (reactants, self['reaction_type'], products)
        
        return result

    def __add__(self,y):
        if isinstance(y,Reaction):
            species = set(self.species()+y.species())
    
            stoic = []
            for spc in species:
                try:
                    stoic.append(dict.__getitem__(self,spc)+y.get(spc, 0))
                except:
                    stoic.append(y[spc]+0)
    
            reaction_type = ('u',self.reaction_type)[self.reaction_type == y.reaction_type]
            
            kwds = dict(zip(species,stoic))
            kwds['reaction_type'] = reaction_type
        elif isinstance(y,Species):
            kwds = self.__add_if_in_spclist(y,self.species())
        else:
            raise TypeError, "Currently, only reactions can be added together"
        
        return Reaction(**kwds)
        
    def __rmul__(self,irrs):
        return self.__mul__(irrs)
        
    def __mul__(self,irrs):
        species = self.species()
        # Using native dict getitem for speed
        stoic = (array([dict.__getitem__(self,k) for k in species], dtype = object)[:,newaxis]*array(irrs).view(ndarray)).swapaxes(0,1)

        # Improvements in Stoic have made this unnecessary
        #stoic = [Stoic(s, s.role) for s in stoic]
        rct_dict = dict(reaction_type = self.reaction_type)
        if not isinstance(irrs, ndarray):
            result = Reaction(**return_updated_dict(dict(zip(species,stoic[0])), rct_dict))
        else:
            result = ReactionArray([Reaction(return_updated_dict(dict(zip(species,stc)), rct_dict)) for stc in stoic], dtype = Reaction)
            
        return result

    def __add_if_in_spclist(self,y,spc_list):
        if not self.has_spc(y):
            raise KeyError, 'Reaction has no components of %s' % y.name
        elif y.name in self.species():
            warn('Already has %s' % y.name)
            return dict(self)
        elif y.exclude:
            raise ValueError, 'Exclude is not supported'
            
        old_species = self.species()

        new_species = [name for name in y.names() if name in spc_list]

        old_stoic = self.stoic()
        new_stoic = sum([self[spc]*y[spc][0] for spc in new_species])

        kwds = dict(zip(old_species+[y.name],old_stoic+[Stoic(new_stoic,'u')]))
        kwds['reaction_type'] = self.reaction_type
            
        return kwds
    
    def get(self, item, default = None):
        try:
            return self.__getitem__(item)
        except:
            return default
            
    def species(self):
        return [k for k, v in self.iteritems() if isinstance(v,Stoic)]
        
    def stoic(self):
        return [self[k] for k in self.species()]
        
    def reactants(self):
        result = [k for k in self.species() if self[k].role == 'r' or (self[k] == 0 and self[k].role == 'u')]
        result += [spc for spc in self.unspecified() if self[spc] < 0]
        
        return result

    def products(self):
        result = [k for k in self.species() if self[k].role == 'p']
        result += [spc for spc in self.unspecified() if self[spc] > 0]
        
        return result

    def unspecified(self):
        return [k for k in self.species() if self[k].role == 'u']

    def has_spc(self,spc_grp):
        return spc_in_list(spc_grp,self.species())
        
    def has_rct(self,spc_grp):
        return spc_in_list(spc_grp,self.reactants())
        
    def has_prd(self,spc_grp):
        return spc_in_list(spc_grp,self.products())
    
    def add_rct_spc(self,y):
        kwds = self.__add_if_in_spclist(y,self.reactants())
        return Reaction(**kwds)
        
    def add_prd_spc(self,y):
        kwds = self.__add_if_in_spclist(y,self.products())
        return Reaction(**kwds)

def spc_in_list(spc_grp,local_list):
    if spc_grp.exclude:
        return not spc_in_list(-spc_grp, local_list)
    else:
        overlapping_spc = set(spc_grp.names()).intersection(local_list)
        return len(overlapping_spc) > 0 

   
class ReactionArray(ndarray):
    """
    ReactionArray is a sub-class of the numpy ndarray that
    supports the Reaction dtype.  This allows indexing
    and mathematical operations on a time-series of Reaction
    objects
    
    Creating an instance follows the typical ndarray __new__ 
    interface.
    
    Only the __getitem__ interface has been implemented; it
    has been implemented to respond appropriately to a Species 
    object.
    """
    def __new__(subtype, *args, **kwds):
        result = array(*args, **kwds).view(type = subtype)
        
        return result

    def __getitem__(self,item):
        if isinstance(item,(Species,str)):
            return array([v[item] for v in self])
        else:
            return ndarray.__getitem__(self,item)

    def get_species_roles(self):
        species_all = self[0].species()
        for hr in self:
            species = hr.reactants()+hr.products()
            if len(species) == len(species_all): break
        else:
            warn('No hour has all species as active.  Ordered species: reactant, unspecified, product.')
            species = [(self[s].sum(), s) for s in hr.species()]
            species.sort()
            species = [s[1] for s in species]
        return species

    def __str__(self):
        species = self.get_species_roles()
        template = '%-21s'+'%11.4f'*(self.shape[0]+1)
        result = ''
        for spc in species:
            result += template % tuple([spc]+self[spc].tolist()+[self[spc].sum()])
            result += '\n'
        
        return result
        

def return_updated_dict(d1, d2):
    d1.update(d2)
    return d1