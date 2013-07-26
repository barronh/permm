import sys
from collections import defaultdict
from warnings import warn
import unittest

from numpy import ndarray, \
                  dtype, \
                  array, \
                  newaxis, \
                  rollaxis, \
                  arange

from PseudoNetCDF.sci_var import PseudoNetCDFVariable

from Species import Species

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
            self.__units[k] = (getattr(kwds[k], 'units', units.get(k, default_unit or 'Unknown'))).strip()
            
    def __getitem__(self, key):
        if isinstance(key, Species):
            values = array([self.data[spc] for spc in key.spc_dict.keys()]).sum(0)
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