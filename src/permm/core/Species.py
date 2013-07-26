from copy import deepcopy
import re
import yaml
from numpy import float32, float64, int8, int16, int32, int64

_spc_def_re = re.compile(r'(?P<stoic>[-+]?[0-9]*\.?[0-9]+)(?P<atom>\S+)(?=\s*\+\s*)?')

def atom_parse(spc_def):
    return dict([(atom, eval(stoic)) for stoic, atom in _spc_def_re.findall(spc_def)])

def atom_guess(spc_name):
    from permm2.mechanisms import atoms as ALL_ATOMS
    if '+' in spc_name or spc_name[:1].isdigit(): return atom_parse(spc_name)
    lastl = ''
    atom_dict = {}
    for l in spc_name:
        if l.isdigit():
            atom_dict[lastl] += int(l) - 1
        elif l in ALL_ATOMS:
            atom_dict.setdefault(l, 0)
            atom_dict[l] += 1
            
        lastl = l
    return atom_dict
    
class Species(object):
    def __init__(self, spc_dict, name = None, exclude = False):
        if isinstance(spc_dict, str):
            if not ':' in spc_dict:
                spc_dict = spc_dict.strip() + ':'
            defs = yaml.load(spc_dict)
            
            if len(defs) > 1:
                raise ValueError('Species class can only initialize one object at a time')
            else:
                (k, v), = defs.items()
                if k is False: k = 'NO'
                name = k
                if v in (None, '', 'GUESS'):
                    atom_dict = atom_guess(k)
                elif v != 'IGNORE':
                    atom_dict = atom_guess(v)
                else:
                    atom_dict = {}
                spc_dict = {k: dict(stoic = 1, atoms = atom_dict)}
        
        self.spc_dict = deepcopy(spc_dict)
        if name:
            self.name = name
        else:
            sep = {False: '+', True: '-'}[exclude]
            prefix = {False: '', True: '-'}[exclude]
            self.name = prefix + sep.join(spc_dict.keys())
            
        for spc, props in self.spc_dict.iteritems():
            props.setdefault('stoic', 1)
            props['role'] = set(props.get('role', 'rup'))
            props.setdefault('atoms', {})

        self.exclude = exclude
    def __getitem__(self, spc_key):
        if isinstance(spc_key, str):
            test_spc = Species({spc_key: dict(stoic = 1, role = set(['r', 'p']), atoms = {})}, name = spc_key)
            return self[test_spc]
        out_spc = {}
        for this_spc, this_props in spc_key.spc_dict.iteritems():
            if this_spc in self.spc_dict:
                check_props = self.spc_dict[this_spc]
                if this_props['role'].issubset(check_props['role']):
                    out_spc[this_spc] = new_props = {}
                    new_props['stoic'] = this_props['stoic'] * check_props['stoic']
                    new_props['atoms'] = deepcopy(this_props['atoms'])
                    new_props['atoms'].update(check_props['atoms'])
                    new_props['role'] = this_props['role']
        if len(out_spc) == 0:
            raise KeyError('%s is not in %s' % (spc_key, self))

        return Species(out_spc, exclude = spc_key.exclude)
    
    def __neg__(self):
        return Species(deepcopy(self.spc_dict), name = '-(%s)' % self.name, exclude = True)
    
    def __str__(self):
        if len(self.spc_dict) == 1:
            (k, v), = self.spc_dict.items()
            if k == self.name and v['stoic'] == 1.:
                return ('%s(%s)' % (self.name, ''.join(self.spc_dict[self.name]['role']))).replace('(pru)', '')
        bool_op = (' = ', ' != ')[self.exclude]
        result = self.name + bool_op + ' + '.join([('%.3f*%s(%s)' % (props['stoic'],spc, ''.join(props['role']))).replace('(pru)', '') for spc, props in self.spc_dict.iteritems()])
        return result
        
    def __repr__(self):
        return self.__str__()
        #try:
        #    exclude = str(self.exclude)
        #except:
        #    exclude = '?'
        #return ndarray.__repr__(self)+'\n      exclude = '+exclude

    def __contains__(self, lhs):
        try:
            self[lhs]
            return True
        except KeyError:
            return False
        
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
        
        if is_number:
            new_props = deepcopy(self.spc_dict)
            for k, props in new_props.iteritems():
                props['stoic'] *= y
            new_name = "%s * %f" % (self.name,float(y))
            new_exclude = y <= 0
            return Species(new_props, name = new_name, exclude = new_exclude)
        else:
            raise TypeError, "Can only multiply species by reactions"

    def __sub__(self, other_species):
        return self.__add__(-other_species)
        
    def __add__(self, other_species):
        return species_sum([self, other_species])

    def atoms(self, atom):
        from permm.mechanisms import atoms as ALL_ATOMS
        if atom in ALL_ATOMS:
            out_props = {}
            for spc, props in self.spc_dict.iteritems():
                mul = props['atoms'].get(atom, 0)
                if mul > 0:
                    out_props[spc] = deepcopy(props)
                    out_props[spc]['stoic'] *= mul
            if out_props == {}:
                raise KeyError, "Atom provided (%s) is not in %s" % (atom, self.name)
            else:
                return Species(out_props, name = self.name + ':' + atom, exclude = self.exclude)
        else:
            raise KeyError, "Atom provided (%s) is not an atom" % atom
    
    def copy(self):
        result =  1*self
        return result

    def stoic(self, spc = None):
        if spc is None:
            spco = self
        else:
            spco = self[spc]
        return sum(v['stoic'] for v in spco.spc_dict.itervalues())

    def iter_species_roles(self):
        for spc, props in self.spc_dict.iteritems():
            for role in props['role']:
                yield spc, role
                
    def role(self, spc = None):
        if spc is None:
            spco = self
        else:
            spco = self[spc]
        roles = [v['role'] for v in spco.spc_dict.itervalues()]

        first_role = roles[0]

        if all([first_role == role for role in roles]):
            return set(first_role)
        else:
            return set('u')
    
    def reactant(self):
        new_props = deepcopy(self.spc_dict)
        for v in new_props.itervalues():
            v['role'] = set('r')
        return Species(new_props, name = self.name, exclude = self.exclude)

    def unspecified(self):
        new_props = deepcopy(self.spc_dict)
        for v in new_props.itervalues():
            v['role'] = set('u')
        return Species(new_props, name = self.name, exclude = self.exclude)

    def product(self):
        new_props = deepcopy(self.spc_dict)
        for v in new_props.itervalues():
            v['role'] = set('p')
        return Species(new_props, name = self.name, exclude = self.exclude)

def species_sum(species_list):
    if not all([isinstance(spc,Species) for spc in species_list]):
        raise TypeError, "Can only add SpeciesGroups"
        
    include_names = set()
    exclude_names = set()
    for next_species in species_list:
        if next_species.exclude:
            exclude_names.update(next_species.spc_dict.keys())
        else:
            include_names.update(next_species.spc_dict.keys())
    
    out_include_names = include_names.difference(exclude_names)
    out_exclude_names = exclude_names.difference(include_names)
    out_props = {}
    if out_include_names != set():
        exclude = False
        out_names = out_include_names
    elif out_exclude_names != set():
        exclude = False
        out_names = out_exclude_names
    else:
        raise ValueError('Sum of species is nothing; check exclusions')
    
    for spc in species_list:
        for spcn, inprops in spc.spc_dict.iteritems():
            if not spcn in out_names: continue

            outprops = out_props.setdefault(spcn, dict(stoic = 0., role = set(), atoms = {}))
            outprops['stoic'] += inprops['stoic']
            outprops['role'].update(inprops['role'])
            outatoms = outprops['atoms']
            inatoms = inprops['atoms']
            for atom, value in inatoms.iteritems():
                try:
                    outatoms[atom] += value
                except:
                    outatoms[atom] = value

    return Species(out_props, exclude = exclude)


import unittest

class SpeciesTestCase(unittest.TestCase):
    def setUp(self):
        self.species = dict(OH = Species(dict(OH = dict(stoic = 1, atoms = dict(H = 1, O = 1))), name = 'OH', exclude = False),
                            HO2 = Species(dict(HO2 = dict(stoic = 1, atoms = dict(H = 1, O = 2))), name = 'HO2', exclude = False),
                            HOx = Species(dict(OH = dict(stoic = 1, atoms = dict(H = 1, O = 1)), HO2 = dict(stoic = 1, atoms = dict(H = 1, O = 2))), name = 'HOx', exclude = False),
                            O3 = Species(dict(O3 = dict(stoic = 1, atoms = dict(O = 3))), name = 'O3', exclude = False),
                            )


    def assertEqualSpecies(self, s1, s2, neg = False):
        self.assertEquals(s1.spc_dict, s2.spc_dict)
        if neg:
            self.assertEquals(s1.exclude, not s2.exclude)
        else:
            self.assertEquals(s1.exclude, s2.exclude)
        
    def testCreate(self):
        s1 = self.species['OH']
        s2 = Species(s1.spc_dict, name = 'hydroxyl')
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
        self.assertEqualSpecies(s3['HO2'], s3[s2])

    def testAtoms(self):
        s2 = self.species['HO2']
        s3 = self.species['HOx']
        self.assertEquals(s2.atoms('O').stoic(s2), 2)
        self.assertEquals(s3.atoms('O').stoic(s2), 2)
        self.assertEquals(s3.atoms('O').stoic(s3), 3)

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
        s4 = Species(s1.spc_dict)
        s5 = s1 + s2 + s3
        self.assertEqualSpecies(s1, s4)
        self.assertEqualSpecies(s1 + s2 + s3, s5)

if __name__ == '__main__':
    unittest.main()
