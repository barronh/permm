from numpy import dtype, \
                  zeros, \
                  ndarray, \
                  array, \
                  float32, float64, int8, int16, int32, int64
import operator

ReactionGroup = str
__all__ = ['Species']

class Species(object):
    """
    Species is a class for species groups.  The simplest case
    is a single species group with a stoichiometry of 1.
    
    The simplest way to create a species is with keywords

        Species(
                name = '<name for the species>',
                names = [<species name1>, ..., <species nameN>],
                stoic = [<factor1>, ..., <factorN>],
                atom_dict = {<species name1> = dict(C=[<factor1>, ..., <factorN>],
                                                N=[<factor1>, ..., <factorN>],
                                                ...)
                         ...
                         <species nameN> = {...}
                         }
                exclude = <boolean>
               )
    
    Species may also be created by combining other species
    with a new name
    
        s1 = Species(**kwds)
        s2 = Species(**kwds)
        Species(s1, s2, name = new_name)
    
    The exclude property determines if the species will be evaluated
    as an exclusion (i.e. not species).  This property may be emulated
    with the __neg__ interface.

    Once a species has been created, it can be 
    1) indexed for stoichiometry
    2) added to other species to create larger species groups,
    3) subtracted from another species
    4) multiplied by a constant
    
    
    1) indexed for stoichiometry
        >>> NO
        Species([(1.0,)], 
              dtype=[('NO', '<f8')])
              exclude = False
        >>> NO['NO']
        array([ 1.])
    
    2) added to other species to create larger species groups,
    2a) inclusive addition
        >>> NO
        Species([(1.0,)], 
              dtype=[('NO', '<f8')])
              exclude = False
        >>> NO2
        Species([(1.0,)], 
              dtype=[('NO2', '<f8')])
              exclude = False
        >>> NOx = NO+NO2
        >>> NOx
        Species([(1.0, 1.0)], 
              dtype=[('NO2', '<f8'), ('NO', '<f8')])
              exclude = False
        
    2b) exclusive addition
        >>> xNO = Species(NO, name='xNO')
        >>> xNO.exclude = True
        >>> xNO
        Species([(1.0,)], 
              dtype=[('NO', '<f8')])
              exclude = True
        >>> NOx
        Species([(1.0, 1.0)], 
              dtype=[('NO2', '<f8'), ('NO', '<f8')])
              exclude = False
        >>> NOx+xNO
        Species([(1.0,)], 
              dtype=[('NO2', '<f8')])
              exclude = False

    2c) exclusive addition with the __neg__ interface
        >>> NOx+-NO
        Species([(1.0,)], 
              dtype=[('NO2', '<f8')])
              exclude = False

    3) subtracted from another species
        >>> NOx
        Species([(1.0, 1.0)], 
              dtype=[('NO2', '<f8'), ('NO', '<f8')])
              exclude = False
        >>> NOx-NO
        Species([(1.0,)], 
              dtype=[('NO2', '<f8')])
              exclude = False
        
    4) multiplied by a constant
        >>> OLE
        Species([(1.0,)], 
              dtype=[('OLE', '<f8')])
              exclude = False
        >>> 2*OLE
        Species([(2.0,)], 
              dtype=[('OLE', '<f8')])
              exclude = False
        >>> Species(2*OLE, name = 'OLEC')
        Species([(2.0,)], 
              dtype=[('OLE', '<f8')])
              exclude = False    
    """
    
    def __init__(self, *args, **kwds):
        if kwds.has_key('names') and kwds.has_key('stoic'):
            # n_names = len(kwds['names'])
            #
            # species_type = dtype(dict(names = kwds['names'], formats = 'd' * n_names))
            #
            # self = zeros((1,), dtype = species_type).view(type = subtype)
            #
            # for name,stoic in zip(kwds['names'],kwds['stoic']):
            #     self[name] += stoic
            self._stoic = {}
            self._stoic.update(dict(zip(kwds['names'], kwds['stoic'])))
            default_atoms = {}
            default_exclude = False
        else:
            args_are_species_groups = all(map(lambda x: isinstance(x, Species),args))

            if args == () or not args_are_species_groups:
                raise ValueError, "Species requires either names and stoic or arguments of type Species"
            else:
                if len(args) == 1:
                    result = args[0]
                else:
                    result = reduce(Species.__add__,args)
            self._stoic = result._stoic.copy()
            self.name = result.name
            default_atoms = result.atom_dict
            default_exclude = result.exclude

        self.atom_dict = kwds.pop('atom_dict', default_atoms).copy()
        try:
            self.name = kwds['name']
        except:
            # Assume that result already has a name
            pass
        self.exclude = kwds.get('exclude', default_exclude)
        self._roles = kwds.get('roles',['u']*len(self.names()))
        
    def __getitem__(self,item):
        """Test"""
        if isinstance(item,Species):
            if item.name in self.keys():
                n = item.name
                return item[n]*self[n]
            else:
                return sum([item[n]*self[n] for n in item.names()])
        else:
            result = dict.__getitem__(self._stoic,item)
            return result
            
    def __str__(self):
        if len(self.names()) == 1:
            if self[self.names()[0]] == 1. and self.name == self.names()[0]:
                return self.name
        bool_op = (' = ', ' != ')[self.exclude]
        result = self.name + bool_op + ' + '.join(['%.3f*%s' % (self[spc],spc) for spc in self.names()])
        return result
        
    def __repr__(self):
        return self.__str__()
        #try:
        #    exclude = str(self.exclude)
        #except:
        #    exclude = '?'
        #return ndarray.__repr__(self)+'\n      exclude = '+exclude
        
    def names(self):
        return [n for n in self.keys()]
    
    def __neg__(self):
        return Species(self, name = '-(%s)' % self.name, exclude = not self.exclude)
    
    def __contains__(self, lhs):
        if isinstance(lhs, Species):
            return len(set(lhs.keys()).intersection(self.keys())) > 0
        elif isinstance(lhs, str):
            return lhs in self.keys()
            
    def __rmul__(self, y):
        return self.__mul__(y)

    def __mul__(self, y):
        is_number = isinstance(y,(int,long,float, float32, float64, int8, int16, int32, int64))
        is_reaction = isinstance(y,ReactionGroup)
        
        if is_number:
            new_name_stoic = [(k, v * y) for k, v in self.iteritems()]
            new_stoic = [v for k, v in new_name_stoic]
            new_name = "%s * %f" % (self.name,float(y))
            new_names = [k for k, v in new_name_stoic]
            new_exclude = y <= 0
            return Species(name = new_name, names = new_names, stoic = new_stoic, exclude = new_exclude, atom_dict = self.atom_dict.copy())
        elif is_reaction:
            raise FutureWarning, "Working on it"
        else:
            raise TypeError, "Can only multiply species by reactions"
    
    def __sub__(self, other_species):
        return self.__add__(-other_species)
        
    def __add__(self, other_species):
        return species_sum([self, other_species])
    
    def atoms(self, atom):
        from mechanisms import atoms
        if atom in atoms:
            name_stoic = [(n, self.atom_dict[n][atom]) for n in self.names() if self.atom_dict.get(n, {}).has_key(atom)]
            if name_stoic == []:
                raise KeyError, "Atom provided (%s) is not in %s" % (atom, self.name)
            else:
                return Species(name = self.name + ':' + atom, names = [n for n, s in name_stoic], stoic = [s for n, s in name_stoic])
        else:
            raise KeyError, "Atom provided (%s) is not an atom" % atom
    
    def copy(self):
        result =  1*self
        result.name = self.name
        return result
    
    def keys(self):
        return self._stoic.keys()

    def iteritems(self):
        return self._stoic.iteritems()

def species_sum(species_list):
    if not all([isinstance(spc,Species) for spc in species_list]):
        raise TypeError, "Can only add SpeciesGroups"
        
    out_species = species_list[0].copy()
    out_names = set(out_species.keys())
    out_name = out_species.name
    out_exclude = out_species.exclude
    out_stoic = {}
    out_atoms = out_species.atom_dict.copy()
    for spc in out_species.keys():
        try:
            out_stoic[spc] += out_species[spc]
        except:
            out_stoic[spc] = float(out_species[spc])
    
    for next_species in species_list[1:]:
        if out_exclude ^ next_species.exclude:
            out_name += ' + ' + next_species.name
            out_names = out_names.difference(next_species.keys())
        else:
            out_name += ' + ' + next_species.name
            out_names.update(next_species.keys())
        
        out_exclude = out_exclude and next_species.exclude
        
        for spc in next_species.keys():
            try:
                out_stoic[spc] += next_species[spc]
            except:
                out_stoic[spc] = next_species[spc]

        if next_species.exclude:
            for spc in next_species.atom_dict.keys():
                del out_atoms[spc]
        else:
            out_atoms.update(next_species.atom_dict)
        
    out_names = list(out_names)
    out_stoic = [out_stoic[name] for name in out_names]
    return Species(name = out_name, names = out_names, \
                    stoic = out_stoic, exclude = out_exclude, atom_dict = out_atoms)

import unittest

class SpeciesTestCase(unittest.TestCase):
    def setUp(self):
        self.species = dict(OH = Species(names = ['OH'], stoic = [1], exclude = False, atom_dict = dict(OH = dict(H = 1, O = 1)), name = 'OH'),
                            HO2 = Species(names = ['HO2'], stoic = [1], exclude = False, atom_dict = dict(HO2 = dict(H = 1, O = 2)), name = 'HO2'),
                            HOx = Species(names = ['HO2', 'OH'], stoic = [1, 1], exclude = False, atom_dict = dict(HO2 = dict(H = 1, O = 2), OH = dict(H = 1, O = 1)), name = 'HOx'),
                            O3 = Species(names = ['O3'], stoic = [1], exclude = False, atom_dict = dict(O2 = dict(O = 3)), name = 'O3'),
                            )


    def assertEqualSpecies(self, s1, s2, neg = False):
        self.assertEquals(s1._stoic, s2._stoic)
        self.assertEquals(s1.atom_dict, s2.atom_dict)
        self.assertEquals(s1._roles, s2._roles)
        if neg:
            self.assertEquals(s1.exclude, not s2.exclude)
        else:
            self.assertEquals(s1.exclude, s2.exclude)
        
    def testCreate(self):
        s1 = self.species['OH']
        s2 = Species(s1, name = 'hydroxyl')
        self.assertEqualSpecies(s1, s2)
        
    def testNeg(self):
        s1 = self.species['OH']
        s2 = -s1
        self.assertEqualSpecies(s1, s2, neg = True)
    
    def testAdd(self):
        s1 = self.species['OH']
        s2 = self.species['HO2']
        s3 = self.species['HOx']
        s4 = s1 + s2
        self.assertEqualSpecies(s3, s4)

    def testSub(self):
        s1 = self.species['OH']
        s2 = self.species['HO2']
        s3 = self.species['HOx']
        s4 = s3 - s2
        self.assertEqualSpecies(s4, s1)

    def testGet(self):
        s2 = self.species['HO2']
        s3 = self.species['HOx']
        self.assertEquals(s3['HO2'], s3[s2])
    
    def testMul(self):
        s1 = self.species['OH']
        s3 = 2 * s1
        s2 = s1 * 2
        self.assertEqualSpecies(s2, s3)
        self.assertRaises(AssertionError, self.assertEqualSpecies, s1, s3)
        
    def testContains(self):
        s1 = self.species['OH']
        s2 = self.species['O3']
        s3 = self.species['HOx']
        self.assertTrue(s1 in s3)
        self.assertFalse(s2 in s3)

    def testCopy(self):
        s1 = self.species['OH']
        s2 = s1.copy()
        self.assertEqualSpecies(s1, s2)

    def testNewFromSpc(self):
        s1 = self.species['OH']
        s2 = self.species['O3']
        s3 = self.species['HOx']
        s4 = Species(s1)
        s5 = Species(s1, s2, s3)
        self.assertEqualSpecies(s1, s4)
        self.assertEqualSpecies(s1 + s2 + s3, s5)

if __name__ == '__main__':
    unittest.main()
