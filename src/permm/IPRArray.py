from collections import defaultdict
from warnings import warn
from numpy import ndarray, \
                  dtype, \
                  array, \
                  newaxis, \
                  rollaxis, \
                  arange
from ProcessGroup import Process
from SpeciesGroup import Species
from PseudoNetCDF.sci_var import PseudoNetCDFVariable
import sys
import unittest

def Processes_ProcDelimSpcDict(proc_names, variables, delim = '_'):
    return dict([(k, Process_ProcDelimSpcDict(k, variables, delim = delim)) for k in proc_names])

def Process_ProcDelimSpcDict(proc_name, variables, delim = '_'):
    """
    Create a single process (proc_name) for all species in the variables
    dictionary where each species/process combination is variable
    with the key <proc_name><delim><spec name>
    proc_name - string process name 
    variables - dictionary with '%s%s%s' % (name, delim, speciesname)
    delim - delimiter
    """
    prefix = '%s%s' % (proc_name, delim)
    lp = len(prefix)
    proc_keys = [k for k in variables.keys() if k[:lp] == prefix]
    process_variables = dict([(k[lp:], variables[k]) for k in proc_keys])
    return Process(proc_name, **process_variables)

class Process(object):
    def __init__(self, name, units = {}, default_unit = None, **kwds):
        self.name = name
        self.names = (name,)
        self.__units = units.copy()
        self.data = {}
        for k in kwds.keys():
            self.data[k] = kwds[k][:]
            self.__units[k] = str.strip(getattr(kwds[k], 'units', units.get(k, default_unit or 'Unknown')))
            
    def __getitem__(self, key):
        if isinstance(key, Species):
            values = array([self.data[spc] for spc in key.names()]).sum(0)
        elif key in self.data.keys():
            values = self.data[key]
        else:
            new_data = dict([(k, v[key]) for k, v in self.data.iteritems()])
            return Process(self.name, self.__units, **new_data)
        return PseudoNetCDFVariable(self, self.name, values.dtype.char, ('TSTEP',), values = values, units = self.__units)

    def __setitem__(self, key, value):
        return self.data.__setitem__(key, value)
    
    def keys(self):
        return self.data.keys()
        
    def set_units(self, k, units):
        if isinstance(k, Species):
            if set([k.name]) == set(k.names()) and not k.exclude:
                self.__units[k.name] = units
            else:
                result = set([self.__units[key] for key in k.names() if self.data.has_key(key)])
                if len(result) == 1:
                    self.__units[k.name] = units
                else:
                    self.__units[k.name] = 'mixed'
        else:
            self.__units[k] = units
        
    def get_units(self, k, defaultunits = None):
        if isinstance(k, Species):
            if self.__units.has_key(k.name):
                return self.__units.get(k.name, defaultunits)
            else:
                result = set([self.__units[key] for key in k.names() if self.data.has_key(key)])
                if len(result) == 1:
                    return result.pop()
                else:
                    return defaultunits or 'mixed'
        else:
            return self.__units[k]
        
    def __getslice__(self, *args, **kwds):
        """
        Use standard ndarray.__getitem__, but retain
        NetRxnArray type
        """
        result = Process(self.name, self.__units, **dict([('%s' % (k), v.__getslice__(*args, **kwds)) for k, v in self.data.iteritems()]))
        return result

    def __generic_math_operator(self, rhs, operation):
        operator = {'+': ndarray.__add__,
                    '-': ndarray.__sub__,
                    '/': ndarray.__div__,
                    '*': ndarray.__mul__
                   }[operation]
        if isinstance(rhs, self.__class__):
            result = Process(name = '%s %s %s' % (self.name, operation, rhs.name))
            result.names = self.names + rhs.names
            for key in self.keys():
                result[key] = self.data[key].copy()
                result.set_units(key, self.get_units(key))
                
            for key in rhs.keys():
                rhs_unit = rhs.get_units(key)
                if result.data.has_key(key):
                    result[key] = operator(result[key], rhs[key])
                    result_unit = result.get_units(key)
                    if rhs_unit != result_unit:
                        result.set_units(key, '%s %s %s' % (result_unit, operation, rhs_unit))
                else:
                    result[key] = rhs[key].copy()
                    result.set_units(key, rhs_unit)
        else:
            result = Process(name = self.name, units = self.__units.copy())
            for key in self.keys():
                try:
                    result[key] = operator(self[key], rhs)
                except:
                    raise TypeError, "It is unclear how to %s IPR and %s; scalars, arrays, species and processes should have meaningful results" % (operation, type(rhs))
            
        return result
        
    def __add__(self,rhs):
        return self.__generic_math_operator(rhs, '+')
        
    def __sub__(self,rhs):
        return self.__generic_math_operator(rhs, '-')
        
    def __div__(self,rhs):
        return self.__generic_math_operator(rhs, '/')
        
    def __mul__(self,rhs):
        return self.__generic_math_operator(rhs, '*')
        
    def sum(self):
        """
        Sum species over time
        """
        result = Process(name = self.name)
        for key in self.keys():
            result[key] = array(self.data[key].sum(), ndmin = 1)
        
        return result
    
class IPR(PseudoNetCDFVariable):
    """
    IPR is a recarray with a field for each species
    and specialized functions for netting, sorting and display
    """
    def __new__(subtype, ipr, spc_names, proc_names, units):
        ipr_dtype = ipr[:].dtype.char
        proc_type = dtype(dict(names = proc_names, formats = ipr_dtype * len(proc_names)))

        spc_type = dtype(dict(names = spc_names, formats = [proc_type] * len(spc_names)))

        result = ipr[:].view(proc_type)[..., 0]

        result = result.view(spc_type)[..., 0]
        
        result = result.view(subtype)
        result.units = units
        result.__dtype = ipr_dtype
        
        return result
        
    def __return_sorted(self, time, comparison_operator, compare_to, reverse = False):
        """
        Internal function for sub-setting species with a binary
        comparison for some time slice
        
        time - slice for the time dimensions
        comparison_operator - a function that evaluates to true 
                              or false when comparing a species
                              netted value to compare_to
        compare_to - a value for comparison by comparison_operator
        """
        
        if isinstance(time, int):
            time = slice(time, time + 1)
        
        # Subset self by time and net for that time
        values = self[time].sum()
        
        # Subset species by comparing netted value to compare_to
        # and create a list of tuples for sorting by value
        species = [(values[k][0], k) for k in self.dtype.names \
                               if comparison_operator(values[k], compare_to)]

        # Sort by value
        species.sort(reverse = reverse)

        # Get names and values for use in rec.fromarrays
        species_names = [k for v, k in species]
        species_values = [v for v, k in species]

        # Construct recarray
        result = rec.fromarrays(species_values, names = species_names)
        
        # Return result
        return result
    def __repr__(self):
        this_dims = str(self.shape)
        this_proc_names = list(set([k for k in self.dtype[0].names]))
        nprocs = len(this_proc_names)
        this_spc_names = list(set([k for k in self.dtype.names]))
        nspcs = len(this_spc_names)
        result = """IPR Array\n  Dimensions: %s\n  Process: %s\n  Species: %s""" % (this_dims, (this_proc_names), str(this_spc_names))
        for prc in this_proc_names:
            result += '\n\n%s:' % prc
            prc_values = self[prc]
            for spc in this_spc_names:
                spc_values = prc_values[spc]
                result += '\n  %s: %s' % (spc, str(spc_values.array()))
        return result
    
    def __str__(self):
        return self.__repr__()

    def __getitem__(self, item):
        """
        Use standard ndarray.__getitem__, but retain
        NetRxnArray type unless item is a species name
        """
        from ProcessGroup import Process
        from SpeciesGroup import Species
        
        if isinstance(item, Process):
            proc_names = set([k for k in item.names])
            this_proc_names = set([k for k in self.dtype[0].names])
            species_names = list(set([k for k in self.dtype.names]))
            
            if item.exclude:
                proc_names = list(this_proc_names.difference(proc_names))
            else:
                proc_names = list(this_proc_names.intersection(proc_names))
            
            spc_type = dtype(dict(names = species_names, formats = self.__dtype * len(species_names)))
    
            result = array([ \
                        array([ndarray.__getitem__(ndarray.__getitem__(self, spc_name), proc_name) for proc_name in proc_names]).sum(0) \
                        for spc_name in species_names
                    ])[newaxis,...]
            result = result.transpose()
            # pre-transpose dims: proc, spc, dim0, dim1, ... dimN
            # post-transpose dims: dimN, ..., dim1, dim0, spc, prc                    

            # Fix N dimension order:
            # for each non-species/process dimension
            #   move the dimension to the front
            # post-fix dims: dim0, dim1, ..., dimN, spc, prc
            for dim in range(1,len(result.shape)-2):
                result = rollaxis(result, dim, start = 0)

            result = IPR(result.copy(), species_names, [item.name], self.units)
            
            
        elif isinstance(item, Species):
            species_names = set([k for k in item.names()])
            this_species_names = set([k for k in self.dtype.names])
            
            if item.exclude:
                species_names = this_species_names.difference(species_names)
            else:
                species_names = this_species_names.intersection(species_names)
            
            if len(species_names) == 0:
                raise KeyError, "%s has no processes" % (item.name,)
                
            proc_size = len(self.dtype[0].names)

            result = array( \
                [ndarray.__getitem__(self, spc_name).copy().view(self.__dtype).reshape(self.shape + (proc_size,))*item[spc_name] for spc_name in species_names]
                ).sum(0)[..., newaxis, :]
                
            species_names = [spc for spc in species_names]
            proc_names = [prc for prc in self.dtype[0].names]
            result = IPR(result.copy(), [item.name], proc_names, units = self.units)
            

        else:
            if isinstance(item, str):
                if item in self.dtype.names:
                    result = self[Species(name = item, names = [item], stoic = [1.], ids = [0])]
                elif item in self.dtype[0].names:
                    result = self[Process(name = item, names = [item])]
                else:
                    raise KeyError, "%s is not a process or a species"

            else:
                result = array(ndarray.__getitem__(self, item), ndmin = 1, copy = True).view(IPR)
                result.__dtype = self.__dtype
                result.units = self.units

        return result

    def __getslice__(self, start, end):
        """
        Use standard ndarray.__getitem__, but retain
        NetRxnArray type
        """
        result = ndarray.__getslice__(self, start, end).view(self.__class__)
        return result

    def __generic_math_operator(self, rhs, operation):
        operator = {'+': ndarray.__add__,
                    '-': ndarray.__sub__,
                    '/': ndarray.__div__,
                    '*': ndarray.__mul__
                   }[operation]
        if isinstance(rhs, self.__class__):
            rhs_dt_sizes = (len(rhs.dtype), len(rhs.dtype[0]))
            self_dt_sizes = (len(self.dtype), len(self.dtype[0]))
            if self_dt_sizes == (1,1):
                if self_dt_sizes == rhs_dt_sizes:
                    spcs = self.dtype.fields.keys()
                    this_prc = self.dtype[0].fields.keys()[0]
                    that_prc = rhs.dtype[0].fields.keys()[0]
                    new_prc_type = dtype(dict(names = ['%s%s%s' % (this_prc, operation, that_prc)], formats = self.__dtype))
                    new_type = dtype(dict(names = spcs, formats = [new_prc_type]*len(spcs)))
                    result = operator(self.array(), rhs.array())
                    
                    result = IPR(result[:, None, None], spc_names = spcs, proc_names = ['%s %s %s' % (this_prc, operation, that_prc)], units = self.units)
                    
                    return result
                else:
                    raise ValueError, "Processes must contain the same number of species and processes"
            else:
                raise ValueError, "Cannot currently %s multi-process processes or species" % (operation,)
        else:
            return operator(self, rhs)

    def __add__(self,rhs):
        return self.__generic_math_operator(rhs, '+')
        
    def __sub__(self,rhs):
        return self.__generic_math_operator(rhs, '-')
        
    def __div__(self,rhs):
        return self.__generic_math_operator(rhs, '/')
        
    def __mul__(self,rhs):
        return self.__generic_math_operator(rhs, '*')
        
    def sum(self):
        """
        Sum net reaction species over time
        """
        
        # Get number of species
        species_size = len(self.dtype.names)
        process_size = len(self.dtype[0].names)

        # Get number of times
        time_size = self.shape[-1]

        # Cast NetRxnArray as ndarray of type float32 and reshape
        # so that species is the last dimension
        result = self.view(self.__dtype).reshape(time_size, species_size, process_size)

        # Sum result across time dimensions using ndarray sum
        # and recast to previous dtype and NetRxnArray type
        result = ndarray.sum(result, 0).view(self.dtype)

        # Return result
        return result
    
    def array(self):
        """
        view IPR array
        """
        result = self.view(type = ndarray, dtype = self.__dtype).view(PseudoNetCDFVariable)
        result.units = self.units
        return result

class TestProcess(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF import PseudoNetCDFFile
        mrg = self.mrg = PseudoNetCDFFile()
        mrg.createDimension('TSTEP', 9)
        mrg.variables = dict(EMIS_NO = PseudoNetCDFVariable(mrg, 'EMIS', 'f', ('TSTEP'), values = arange(9), units = 'ppb'),
                             CHEM_NO = PseudoNetCDFVariable(mrg, 'CHEM', 'f', ('TSTEP'), values = arange(9)),
                             EMIS_NO2 = PseudoNetCDFVariable(mrg, 'EMIS', 'f', ('TSTEP'), values = arange(9), units = 'ppb'),
                             CHEM_NO2 = PseudoNetCDFVariable(mrg, 'CHEM', 'f', ('TSTEP'), values = arange(9), units = 'ppb')
                            )
        self.processes = {}
        self.species = dict(NO = Species(name = 'NO', stoic = [1], names = ['NO']),
                            NO2 = Species(name = 'NO2', stoic = [1], names = ['NO2']),
                            )
        exec('NOx = NO + NO2', None, self.species)
        self.testProcessFromCMAQ()

    def testProcessFromCMAQ(self):
        self.processes.update(Processes_ProcDelimSpcDict('EMIS CHEM'.split(), self.mrg.variables))
        self.assertTrue(set(self.processes.keys()) == set(['EMIS', 'CHEM']))
    
    def testIndex(self):
        testproc = self.processes['EMIS']
        NO = self.species['NO']
        NO2 = self.species['NO2']
        NOx = self.species['NOx']
        self.assertTrue((testproc[NO] == arange(9)).all())
        self.assertTrue((testproc[NOx] == (2 * arange(9))).all())
        self.assertTrue((testproc[:3][NO] == arange(9)[:3]).all())
        self.assertTrue((testproc[NO][:3] == arange(9)[:3]).all())

    def testAdd(self):
        NO = self.species['NO']
        NO2 = self.species['NO2']
        NOx = self.species['NOx']
        testproc1 = self.processes['EMIS']
        testproc2 = self.processes['CHEM']
        testproc3 = testproc1 + testproc2
        self.assertTrue((testproc3[NO] == (arange(9) * 2)).all())
        self.assertTrue((testproc3[NO2] == (arange(9) * 2)).all())
        self.assertTrue((testproc3[NOx] == (arange(9) * 4)).all())

    def testSum(self):
        NO = self.species['NO']
        NO2 = self.species['NO2']
        NOx = self.species['NOx']
        testproc1 = self.processes['EMIS']
        testproc2 = self.processes['CHEM']
        testproc3 = testproc1 + testproc2
        self.assertTrue(testproc2.sum()[NO] == arange(9).sum())
        self.assertTrue(testproc2[NO].sum() == arange(9).sum())
        self.assertTrue((testproc3.sum()[NO2] == (2 * arange(9).sum())))
        self.assertTrue(testproc3[NOx].sum() == (arange(9) * 4).sum())
    
    def testUnits(self):
        NO = self.species['NO']
        self.assertTrue(self.processes['EMIS'].get_units(NO) == 'ppb')
        self.assertTrue(self.processes['CHEM'].get_units(NO) == 'Unknown')

if __name__ == '__main__':
    unittest.main()