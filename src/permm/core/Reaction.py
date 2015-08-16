import operator
import re
from warnings import warn
from copy import deepcopy
from numpy import ndarray, \
                  newaxis, \
                  float64, \
                  vectorize, \
                  rollaxis, \
                  array, \
                  indices, \
                  equal, \
                  greater, \
                  less
from numpy.ma import sum, masked_less, masked_greater

from permm.core.Species import Species

ReactionGroup = str

__all__ = ['Stoic', 'Reaction']


class Stoic(ndarray):
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
    __array_priority__ = 1001.
    def __new__(subtype, f, role):
        result = array(f).view(subtype)
        result.role = role
        return result
    
    def __array_wrap__(self, obj, context = None):
        result = obj.view(Stoic)
        if context is not None:
            args = context[1]
            if len(args) > 1:
                lhs, rhs = args[:2]
                if hasattr(lhs, 'role'):
                    if hasattr(rhs, 'role'):
                        if lhs.role != rhs.role:
                            result.role = 'u'
                        else:
                            result.role = lhs.role
                    else:
                            result.role = lhs.role                        
                else:
                    result.role = 'u'

        return result
        
    def __repr__(self):
        return "Stoic(%s, role = '%s')" % (str(self.asarray()), getattr(self, 'role', 'None'))
    
    def __str__(self):
        return self.__repr__()
            
    def asarray(self):
        return self.view(ndarray)
        

StoicArray = vectorize(Stoic, otypes = [object])
def AddRole(s1, s2):
    if isinstance(s1,Stoic):
        return s1
    else:
        return Stoic(s1, s2.role)
        
def StoicAdd(s1, s2):
    return Stoic.__add__(AddRole(s1,s2), s2)

StoicAdd = vectorize(StoicAdd, otypes = [object])

def ParseReactionString(rxn_str):
    """
    ReactionFromString is a convenience function.  It creates
    Reaction objects from a string with the following pattern:
        "(?P<reactants>.*)=(?P<rxn_type>[kj])[>]\s*(?P<products>.*)"
    where each reactant or product matches the following pattern:
        "(\s?(?P<sign>[+-])?\s?)?((?P<stoic>\d(\.(\d{1,3}(E\d{2})?)?)?)\*)?(?P<name>[A-Z]\w{0,3})\s?"

    For example:
        OH + OLE =k> 0.8*FORM + 0.33*ALD2 + 0.62*ALDX + 0.8*XO2 + 0.95*HO2 - 0.7 PAR
    """
    stoics = {}
    
    reaction_re = re.compile("(?P<reactants>.*)(?:=|->\[(?P<rxn_type>[kjdeu])\])\s*(?P<products>.*)")

    species_re = re.compile("(\s?(?P<sign>[+-])?\s?)?((?P<stoic>\d{0,1}(\.(\d{1,3}(E\d{2})?)?)?)\*)?(?P<name>[a-zA-Z]\w*)(?:[ +=]|$)+",re.M)
    
    reaction_match = reaction_re.match(rxn_str)
    if reaction_match is None:
        raise SyntaxError, "Reactions must match the following patter\n\n<%%d> stoic*spc ([+-] stoic*spc)* =[kjude]> [stoic*spc] ([+-] stoic*spc)*\n\n%s" % (rxn_str,)

    reactants = reaction_match.groupdict()['reactants']
    reaction_type = reaction_match.groupdict()['rxn_type']
    if reaction_type is None:
        reaction_type = 'u'
    products = reaction_match.groupdict()['products']

    for spc in species_re.finditer(reactants):
        spc_g = spc.groupdict()
        name = spc_g['name']
        sign = spc_g['sign']
        stoic = spc_g['stoic']
        if sign is None:
            sign = ''

        if stoic is None:
            stoic = '1'

        if (name, 'r') in stoics:
            stoics[name, 'r'] += -float(sign + stoic)
        else:
            stoics[name, 'r'] = -float(sign + stoic)

    
    for spc in species_re.finditer(products):
        spc_g = spc.groupdict()
        name = spc_g['name']
        sign = spc_g['sign']
        stoic = spc_g['stoic']
        if sign is None:
            sign = '+'

        if stoic is None:
            stoic = '1'
        
        value = float(sign + stoic)
        if (name, 'p') in stoics:
            stoics[name, 'p'] += value
        else:
            stoics[name, 'p'] = value

    return stoics, reaction_type
    
def ReactionFromString(rxn_str):
    return Reaction(*ParseReactionString(rxn_str))

class Reaction(object):
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
    __array_priority__ = 1000.

    def __init__(self, stoic, reaction_type = 'k', safe = False):
        """
            roles - dictionary of specific roles; for species whose 
                    stoichiometric sign is inconsistent with their 
                    role (i.e. X + Y =k> Z - .6*PAR)
            reaction_type - 'k' is kinetic, 'j' is photolysis, 'n' is net
            stoic - stoichiometry values provided as keywords by species; values
                    can be scalars or ndarrays
            safe - when species roles are not found, find any version of
                   that species to prevent an error
        """
        if isinstance(stoic, str):
            stoic, reaction_type = ParseReactionString(stoic)

        self.reaction_type = reaction_type
        self._safe = safe
        self._stoic = deepcopy(stoic)
        self._update_roles()
        try:
            first_key = self._stoic.keys()[0]
            self.shape = self._stoic[first_key].shape
        except AttributeError, (e):
            self.shape = ()
            for k, v in self._stoic.iteritems():
                self._stoic[k] = float64(v)
        except IndexError, (e):
            self.shape = ()
            
    
    def _update_roles(self):
        keys = self._stoic.keys()
        self._species = tuple(set([spcn for spcn, role in keys]))
        self._reactants = tuple(set([spcn for spcn, role in keys if role == 'r']))
        self._products = tuple(set([spcn for spcn, role in keys if role == 'p']))
        self._unspecified = tuple(set([spcn for spcn, role in keys if role == 'u']))
        
    def roles(self, spc):
        roles = set()
        if isinstance(spc, str):
            if spc in self._reactants:
                roles.update('r')
            if spc in self._products:
                roles.update('p')
            if spc in self._unspecified:
                roles.update('u')
                
        elif isinstance(spc, Species):
            roles.update(*[self.roles(spcn) for spcn in spc.spc_dict.iterkeys()])
        
        return roles
            
    def __contains__(self, lhs):
        """
        Test if reaction has a species
        """
        if isinstance(lhs, Species):
            return self.has_spc(lhs)
        elif isinstance(lhs, str):
            return lhs in self._species
        else:
            raise TypeError, 'Unknown comparison: __contains__ for Reaction and %s' % str(type(lhs))
            
    def __getitem__(self, item):
        """
        Return a stoichiometry for a species or, if item is not a species, return
        a new reaction where stoichiometry is subset using item
        """
        if isinstance(item, Species):
            return self.get_spc(item)
        elif isinstance(item, str):
            if  item in self._species:
                outvalue = 0.
                outrole = set()
                for role in self.roles(item):
                    outvalue += self._stoic[item, role]
                    outrole.update(role)

                if len(outrole) == 1:
                    role, = outrole
                else:
                    role = 'u'
                
            return Stoic(outvalue, role = role)
        else:
            return Reaction(dict([(k,v[item]) for k, v in self._stoic.iteritems()]), reaction_type = self.reaction_type)
    
    def __str__(self):
        """
        Report all values followed by the sum value for the reaction
        """
        result = ''
        temp = self.copy()
        temp.shape = ()
        if self.shape != ():
            result += '%d Reactions with shape %s\n' % (reduce(int.__mul__, self.shape), str(self.shape))
            for idx in indices(self.shape).reshape(len(self.shape), -1).swapaxes(0,1):
                idx = tuple(idx)
                result += str(idx) + ': '
                for spc in self._stoic.keys():
                    temp._stoic[spc] = self._stoic[spc][idx]
                temp._update_roles()
                result += str(temp)+', \n'
            result = result[:-3]+'\n'
        sum_result = self.display(digits = 5, nspc = 1000)

        if result != '':
            result += '-' * (len(sum_result)+5) + '\nSum: '
        result += sum_result
            
        return result

    def __repr__(self):
        """
        Representation is same as string
        """
        return self.__str__()

    def display(self, digits = 5, nspc = 3, formatter = 'g'):
        reactants = []
        products = []
        for (spcn, role), v in self._stoic.iteritems():
            if role == 'r':
                reactants.append((self._stoic[spcn, role].sum(), spcn))
            elif role == 'p':
                products.append((self._stoic[spcn, role].sum(), spcn))
            else:
                v = self._stoic[spcn, role].sum()
                if v < 0:
                    reactants.append((v, spcn))
                else:
                    products.append((v, spcn))

        reactants.sort(reverse=False)
        products.sort(reverse=True)

        if digits == None:
            str_temp = '%s'
            reactant_str = ' + '.join([str_temp % rct for stoic, rct in reactants][:nspc])
            product_str = ' + '.join([str_temp % prd for stoic, prd in products][:nspc])
        elif digits == -1:
            reactant_str = ' + '.join(['*'.join([str(-stoic), str(rct)]) for stoic, rct in reactants][:nspc])
            product_str = ' + '.join(['*'.join([str(stoic), str(prd)]) for stoic, prd in products][:nspc])
        else:
            str_temp = '%%.%d%s*%%s' % (digits, formatter)
            reactant_str = ' + '.join([str_temp % (-stoic,rct) for stoic, rct in reactants][:nspc])
            product_str = ' + '.join([str_temp % (stoic,prd) for stoic, prd in products][:nspc])

        if len(reactants) > nspc:
            reactant_str += ' + ...'
        if len(products) > nspc:
            product_str += ' + ...'

        return '%s ->[%s] %s' % (reactant_str, self.reaction_type, product_str)
        
    def __add__(self,rhs):
        """
        Add reactions to make a net reaction or add species to an existing reaction.
        """
        if isinstance(rhs,Reaction):
            kwds = {}
            
            for key, v in self._stoic.iteritems():
                kwds[key] = self._stoic[key].copy()
            
            for key, v in rhs._stoic.iteritems():
                if kwds.has_key(key):
                    kwds[key] += rhs._stoic[key].copy()
                else:
                    kwds[key] = rhs._stoic[key].copy()

            if self.reaction_type == rhs.reaction_type:
                reaction_type = self.reaction_type
            else:
                reaction_type = 'n'

        elif isinstance(rhs,Species):
            return self.__add_if_in_spclist(rhs,self._species)

        else:
            raise TypeError, "Currently, only reactions can be added together"
        
        return Reaction(kwds, reaction_type = reaction_type)
        
    def __rmul__(self,irrs):
        return self.__mul__(irrs)
        
    def __mul__(self,irrs):
        species = self._species
        values = dict([(k,v*irrs) for k, v in self._stoic.iteritems()])
        result = Reaction(values, reaction_type = self.reaction_type)
        return result

    def __add_if_in_spclist(self,rhs,spc_list):
        if not self.has_spc(rhs):
            raise KeyError, 'Reaction has no components of %s' % rhs.name
        elif rhs.name in self._species:
            warn('Already has %s' % rhs.name)
            return self.copy()
        elif rhs.exclude:
            raise ValueError, 'Exclude is not supported'
        result = self.copy()
        new_stoic = result[rhs]
        
        result._roles[rhs.name] = new_stoic.role
        result._species += (rhs.name,)
        result._stoic[rhs.name] = new_stoic.view(ndarray)
        result._update_roles()

        return result
            
    def copy(self):
        """
        Create a copy of the reaction such that stoichiometry 
        are not shared
        """
        return Reaction(dict([(k, v.copy()) for k, v in self._stoic.iteritems()]), reaction_type = self.reaction_type, )

    def sum(self, axis = None):
        """
        Sum stoichiometries and create a scalar reaction
        """
        result = self.copy()
        for key in result._stoic.keys():
            result._stoic[key] = self._stoic[key].sum(axis)
        result.shape = ()
        result._update_roles()
        return result
        
    def mean(self, axis = None):
        """
        Mean stoichiometries and create a scalar reaction
        """
        result = self.copy()
        for spc in result._species:
            result._stoic[spc] = self._stoic[spc].mean(axis)
        result.shape = ()
        result._update_roles()
        return result
        
    def reactants(self):
        """
        Report all species acting as reactants; includes negative
        stoichiometries for unspecified role species
        """
        return self._reactants
        
    def products(self):
        """
        Report all species acting as products; includes positive
        stoichiometries for unspecified role species
        """
        return self._products
        
    def species(self):
        """
        Report all species
        """
        return self._species

    def unspecified(self):
        """
        Report all species with unspecified roles
        """
        return self._unspecified

    def get(self, item, default = None):
        try:
            return self.__getitem__(item)
        except:
            return default

    def get_spc(self, *args):
        nargs = len(args)
        if nargs > 2:
            raise TypeError('get_spc expected at most 2 arguments, got %d' % nargs)
            
        item = args[0]

        values = []
        roles = []
        if item.exclude:
            item_spc_roles = [k for k in item.iter_species_roles()]
            item_spcs = [k[0] for k in item_spc_roles]
            for spc, role in self._stoic:
                if role == 'u' and spc in item_spcs and (spc, 'u') not in item_spc_roles:
                    val = self._stoic[spc, role]
                    if greater(val, 0).all():
                        role = 'p'
                    elif less(val, 0).all():
                        role = 'r'
                
                if not (spc, role) in item_spc_roles:
                    values.append(1 * self._stoic[spc, role])
                    roles.append(role)
                    
                
        else:
            for spc, props in item.spc_dict.iteritems():
                spc_roles = props['role']
                for role in spc_roles:
                    if (spc, role) in self._stoic:
                        values.append(item.spc_dict[spc]['stoic'] * self._stoic[spc, role])
                        roles.append(role)
                
                if 'u' not in spc_roles and (spc, 'u') in self._stoic:
                    val = self._stoic[spc, 'u']
                    if role == 'p':
                        values.append(item.spc_dict[spc]['stoic'] * masked_less(val, 0).filled(0))
                        roles.append(role)
                    elif role == 'r':
                        values.append(item.spc_dict[spc]['stoic'] * masked_greater(val, 0).filled(0))
                        roles.append(role)
                    
                
                
        if len(values) == 0:
            if self._safe:
                if nargs == 1:
                    if item.name in self:
                        warn("%s not in %s; trying all roles" % (item, self.sum()))
                        return self.get_spc(Species(item.name, exclude = item.exclude))
                    raise KeyError, "%s does not contain %s" % (str(self.sum()), str(item))
                else:
                    return args[1]
            else:
                raise KeyError, "%s does not contain %s" % (str(self.sum()), str(item))
            
                
        last_role = roles[-1]
        same_role = all([last_role == role for role in roles])
        
        if same_role:
            role = last_role
        else:
            role = 'u'
        return Stoic(sum(values, axis = 0), role = role)
    
    def produces(self, item):
        try:
            return self.get_spc(item.product())
        except KeyError:
            return 0 * self._stoic[self._stoic.keys()[0]]

    def consumes(self, item):
        try:
            return -self.get_spc(item.reactant())
        except KeyError:
            return 0 * self._stoic[self._stoic.keys()[0]]
    
    def has_spc(self,spc_grp):
        for spc_role in spc_grp.iter_species_roles():
            if spc_role in self._stoic:
                return True != spc_grp.exclude
        return False != spc_grp.exclude

    def has_rct(self,spc_grp):
        for name in spc_grp.names():
            if not spc_grp.contains_species_role(name, 'r'):
                raise TypeError('Requesting reactant role from species %s with %s that has roles %s' % (spc_grp.name, name, str(list(spc_grp.spc_dict[name]['role']))))
        return self.has_spc(spc_grp.reactant())
        
    def has_prd(self,spc_grp):
        for name in spc_grp.names():
            if not spc_grp.contains_species_role(name, 'p'):
                raise TypeError('Requesting product role from species %s with %s that has roles %s' % (spc_grp.name, name, str(list(spc_grp.spc_dict[name]['role']))))
        return self.has_spc(spc_grp.product())
    
    def has_unspc(self,spc_grp):
        for name in spc_grp.names():
            if not spc_grp.contains_species_role(name, 'u'):
                raise TypeError('Requesting unspecified role from species %s with %s that has roles %s' % (spc_grp.name, name, str(list(spc_grp.spc_dict[name]['role']))))
        return self.has_spc(spc_grp.unspecified())
    
    def net(self, spco = None):
        if spco is None:
            spcs_to_condense = self._species
        else:
            if spco.name in self:
                spcs_to_condense = [spco.name]
            else:
                spcs_to_condense = spco.spc_dict.keys()
        net_spcs = []
        for spc in spcs_to_condense:
            prd = spc in self._products
            rct = spc in self._reactants
            unk = spc in self._unspecified
            if not prd ^ rct ^ unk:
                net_spcs.append(spc)
        kwds = deepcopy(self._stoic)
        dummy_stoic = Stoic(0., role = 'u')
        for spc in net_spcs:
            rval = kwds.pop((spc, 'r'), dummy_stoic)
            pval = kwds.pop((spc, 'p'), dummy_stoic)
            uval = kwds.pop((spc, 'u'), dummy_stoic)
            new_val = rval + pval + uval

            if not equal(new_val, 0.).all():
                kwds[spc, 'u'] = new_val

        return Reaction(kwds, reaction_type = self.reaction_type)

    def condense(self, spco, name = None):
        if not isinstance(spco, Species):
            raise TypeError('condense can only take Species objects')
        if spco in self:
            net_keys = [(spc, role) for spc, role in spco.iter_species_roles() if (spc, role) in self._stoic]
            new_roles = set([role for spc, role in net_keys])
            new_val_dict = dict([(role, []) for role in new_roles])
            new_name = name or spco.name
            
            kwds = deepcopy(self._stoic)
            for spc, role in net_keys:
                new_val_dict[role].append(kwds.pop((spc, role)))

            for new_role, new_vals in new_val_dict.iteritems():
                kwds[new_name, new_role] = reduce(operator.add, new_vals)

            return Reaction(kwds, reaction_type = self.reaction_type)
        else:
            return self.copy()
        
import unittest

class ReactionTestCase(unittest.TestCase):
    def setUp(self):
        self.spcs = dict(NO2 = Species('NO2'),
                         NO = Species('"NO"'),
                         HNO3 = Species('HNO3'),
                         O = Species('O'),
                         OH = Species('OH'),
                         HO2 = Species('HO2'),
                         O3 = Species('O3'),
                         PAR = Species('PAR: C'),
                         NTR = Species('NTR: IGNORE'),
                         FORM = Species('FORM: CH2O'),
                         ALD2 = Species('ALD2: CH3CHO'),
                         ALDX = Species('ALDX: CH3CH2CHO'),
                         )
        exec('NOx = NO2 + NO', globals(), self.spcs)
        exec('ALD = ALDX + ALD2 + FORM', globals(), self.spcs)
        rxn_strings = {'NO2hv': 'NO2 ->[j] NO + O',
                       'OplO3': 'O + O2 + M ->[k] O3 + M',
                       'NTRplOH': 'NTR + OH ->[k] HNO3 + HO2 + 0.330*FORM + 0.330*ALD2 + 0.330*ALDX - 0.660*PAR'
                      }
        self.rxns = {}
        for label, rxn_str in rxn_strings.iteritems():
            self.rxns[label] = Reaction(rxn_str)
        
    def testHasSpc(self):
        r1 = self.rxns['NO2hv']
        spcs = [self.spcs[k] for k in ['NO2', 'NO', 'NOx']]
        for spc in spcs:
            self.assertTrue(r1.has_spc(spc))
            
    def testHasRole(self):
        r1 = self.rxns['NO2hv']
        rs = self.spcs['NO2']
        ps = self.spcs['NO']
        self.assertTrue(r1.has_rct(rs))
        self.assertFalse(r1.has_rct(ps))
        self.assertFalse(r1.has_prd(rs))
        self.assertTrue(r1.has_prd(ps))
        
    def testConsumes(self):
        r1 = self.rxns['NO2hv']
        rs = self.spcs['NO2']
        ps = self.spcs['NO']
        ns = self.spcs['O3']
        self.assertTrue(r1.consumes(rs) == 1)
        self.assertTrue(r1.consumes(ps) == 0)
        self.assertTrue(r1.consumes(ns) == 0)
        
    def testProduces(self):
        r1 = self.rxns['NO2hv']
        rs = self.spcs['NO2']
        ps = self.spcs['NO']
        ns = self.spcs['O3']
        self.assertTrue(r1.produces(rs) == 0)
        self.assertTrue(r1.produces(ps) == 1)
        self.assertTrue(r1.produces(ns) == 0)
        
    def testGet(self):
        r1 = self.rxns['NO2hv']
        rs = self.spcs['NO2']
        self.assertTrue(r1[rs] == -1)
        
    def testReactants(self):
        r1 = self.rxns['NTRplOH']
        self.assertEquals(set(r1.reactants()), set(('OH', 'NTR')))
        
        
    def testProducts(self):
        r1 = self.rxns['NTRplOH']
        self.assertEquals(set(r1.products()), set(('PAR', 'FORM', 'ALD2', 'HNO3', 'HO2', 'ALDX')))
        
    def testMulScalar(self):
        r1 = self.rxns['NTRplOH']
        a = 5.83
        r2 = r1 * a

        for spcn in "NTR OH".split():
            spc = self.spcs[spcn]
            self.assertAlmostEqual(r2[spc], -a)

        for spcn in "HNO3 HO2".split():
            spc = self.spcs[spcn]
            self.assertAlmostEqual(r2[spc], a)

        for spcn in "FORM ALD2 ALDX".split():
            spc = self.spcs[spcn]
            self.assertAlmostEqual(r2[spc], a * .33)

        for spcn in ["ALD"]:
            spc = self.spcs[spcn]
            self.assertAlmostEqual(r2[spc], a * .99)

        for spcn in ["PAR"]:
            spc = self.spcs[spcn]
            self.assertAlmostEqual(r2[spc], -a * .66)
        
    def testMulArray(self):
        from numpy import arange, round
        a = arange(0, 60, dtype = 'd').reshape(3,4,5) + .3
        r1 = self.rxns['NTRplOH']
        r2 = r1 * a

        for spcn in "NTR OH".split():
            spc = self.spcs[spcn]
            self.assertTrue((r2[spc] == -a).all())

        for spcn in "HNO3 HO2".split():
            spc = self.spcs[spcn]
            self.assertTrue((r2[spc] == a).all())

        for spcn in "FORM ALD2 ALDX".split():
            spc = self.spcs[spcn]
            self.assertTrue((r2[spc] == (a * float64(.33))).all())

        for spcn in ["ALD"]:
            spc = self.spcs[spcn]
            # Value is not exact because of addition
            self.assertTrue((round(r2[spc], decimals = 8) == round(a * float64(.99), decimals = 8)).all())

        for spcn in ["PAR"]:
            spc = self.spcs[spcn]
            self.assertTrue((r2[spc] == (-a * .66)).all())
        
if __name__ == '__main__':
    unittest.main()
    r1 = Reaction('NO2 ->[j] NO + O')
    r2 = Reaction('O + O2 + M ->[k] O3 + M')
    NO2 = Species('NO2')
    NO = Species('NO')
    O = Species('O')
    O3 = Species('O3')
    Ox = O3 + NO2
    Ox.name = 'Ox'
    print r1[NO2]
    print r1[NO]
    print r1[O]
    from numpy import arange
    nr = (r1 * arange(9) * 1.1 + r2 * arange(9)[::-1]).net()
    print nr[:4][O.reactant()]
    print nr[O]
    try:
        print nr[O.reactant()]
    except:
        pass
    print nr.net()
    print nr.net().condense(Ox)