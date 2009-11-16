from numpy import dtype, \
                  zeros, \
                  ndarray, \
                  array

import operator

ReactionGroup = str
__all__ = ['Species']

class Species(ndarray):
    """
    Species is a class for species groups.  The simplest case
    is a single species group with a stoichiometry of 1.
    
    The simplest way to create a species is with keywords

        Species(
                name = '<name for the species>',
                names = [<species name1>, ..., <species nameN>],
                stoic = [<factor1>, ..., <factorN>],
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
    
    def __new__(subtype, *args, **kwds):
        if kwds.has_key('names') and kwds.has_key('stoic'):
            n_names = len(kwds['names'])

            species_type = dtype(dict(names = kwds['names'], formats = 'd' * n_names))
            
            result = zeros((1,), dtype = species_type).view(type = subtype)

            for name,stoic in zip(kwds['names'],kwds['stoic']):
                result[name] += stoic
                
        else:
            args_are_species_groups = reduce(operator.and_,[True]+map(lambda x: isinstance(x, Species),args))

            if args == () or not args_are_species_groups:
                raise ValueError, "Species requires either names and stoic or arguments of type Species"
            else:
                if len(args) == 1:
                    arg = args[0]
                    result = Species(name = arg.name, names = arg.names(), stoic = [arg[n] for n in arg.names()])
                else:
                    result = reduce(operator.add,args)

        result.name = kwds['name']
        result.exclude = kwds.get('exclude', False)
        result.__roles = kwds.get('roles',['u']*len(result.names()))
        return result
    
    def __getitem__(self,item):
        if isinstance(item,Species):
            if self.has_key(item.name):
                return item[n]*self[n]
            else:
                return sum([item[n]*self[n] for n in item.names()])
        else:
            result = ndarray.__getitem__(self,item)
            return result.view(ndarray)
            
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
        return [n for n in self.dtype.names]
    
    def __neg__(self):
        return Species(self, name = 'x%s' % self.name, exclude = not self.exclude)
    
    def __contains__(self, lhs):
        if isinstance(lhs, Species):
            return len(set(lhs.dtype.names).intersection(self.dtype.names)) > 0
        elif isinstance(lhs, str):
            return lhs in self.dtype.names
            
    def __rmul__(self, y):
        return self.__mul__(y)

    def __mul__(self, y):
        is_number = isinstance(y,(int,long,float))
        is_reaction = isinstance(y,ReactionGroup)
        
        if is_number:
            new_stoic = array(self.tolist()[0]) * y
            new_name = "%s * %f" % (self.name,float(y))
            new_names = [n for n in self.names()]
            new_exclude = y <= 0
            return Species(name = new_name, names = new_names, stoic = new_stoic, exclude = new_exclude)
        elif is_reaction:
            raise FutureWarning, "Working on it"
        else:
            raise TypeError, "Can only multiply species by reactions"
    
    def __sub__(self, other_species):
        return self.__add__(-other_species)
        
    def __add__(self, other_species):
        return species_sum([self, other_species])

def species_sum(species_list):
    if not all([isinstance(spc,Species) for spc in species_list]):
        raise TypeError, "Can only add SpeciesGroups"
        
    out_species = species_list[0]
    out_names = set(out_species.dtype.names)
    out_name = out_species.name
    out_exclude = out_species.exclude
    out_stoic = {}
    for spc in out_species.dtype.names:
        try:
            out_stoic[spc] += out_species[spc]
        except:
            out_stoic[spc] = out_species[spc].copy()

    for next_species in species_list[1:]:
        if out_exclude ^ next_species.exclude:
            out_name += ' + ' + next_species.name
            out_names = out_names.difference(next_species.dtype.names)
        else:
            out_name += ' + ' + next_species.name
            out_names.update(next_species.dtype.names)
        
        out_exclude = out_exclude and next_species.exclude
        
        for spc in next_species.dtype.names:
            try:
                out_stoic[spc] += next_species[spc]
            except:
                out_stoic[spc] = next_species[spc]
    
    out_names = list(out_names)
    out_stoic = [out_stoic[name] for name in out_names]
    return Species(name = out_name, names = out_names, \
                    stoic = out_stoic, exclude = out_exclude)
