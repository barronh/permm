from numpy import ndarray, \
                  dtype, \
                  array, \
                  newaxis
from ProcessGroup import Process
from SpeciesGroup import Species
import sys

class IPR(ndarray):
    """
    IPR is a recarray with a field for each species
    and specialized functions for netting, sorting and display
    """
    def __new__(subtype, ipr, spc_names, proc_names):
        proc_type = dtype(dict(names = proc_names, formats = 'f' * len(proc_names)))

        spc_type = dtype(dict(names = spc_names, formats = [proc_type] * len(spc_names)))

        result = ipr.view(proc_type)[:, :, 0]

        result = result.view(spc_type)[:, 0]
        
        result = result.view(subtype)
        
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
            
            spc_type = dtype(dict(names = species_names, formats = 'f' * len(species_names)))
    
            result = array([ \
                        array([ndarray.__getitem__(ndarray.__getitem__(self, spc_name), proc_name) for proc_name in proc_names]).sum(0) \
                        for spc_name in species_names
                    ], ndmin = 3).swapaxes(0, 2)
            result = IPR(result.copy(), species_names, [item.name])
            
            
        elif isinstance(item, Species):
            species_names = set([k for k in item.names()])
            this_species_names = set([k for k in self.dtype.names])
            
            if item.exclude:
                species_names = this_species_names.difference(species_names)
            else:
                species_names = this_species_names.intersection(species_names)
            
            if len(species_names) == 0:
                raise KeyError, "%s has no processes" % (item.name,)
                
            time_size = self.shape[0]
            proc_size = len(self.dtype[0].names)

            result = array( \
                [ndarray.__getitem__(self, spc_name).copy().view('f').reshape(time_size, proc_size)*item[spc_name][0] for spc_name in species_names]
                ).sum(0)[:, newaxis, :]
                
            species_names = [spc for spc in species_names]
            proc_names = [prc for prc in self.dtype[0].names]
            result = IPR(result.copy(), [item.name], proc_names)
            

        else:
            if isinstance(item, str):
                if item in self.dtype.names:
                    result = self[Species(name = item, names = [item], stoic = [1.], ids = [0])]
                elif item in self.dtype[0].names:
                    result = self[Process(name = item, names = [item])]
                else:
                    raise KeyError, "%s is not a process or a species"

            else:
                result = ndarray.__getitem__(self, item)
                result = result.view(self.__class__)

        return result

    def __getslice__(self, start, end):
        """
        Use standard ndarray.__getitem__, but retain
        NetRxnArray type
        """
        result = ndarray.__getslice__(self, start, end).view(self.__class__)
        return result

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
        result = self.view('f').reshape(time_size, species_size, process_size)

        # Sum result across time dimensions using ndarray sum
        # and recast to previous dtype and NetRxnArray type
        result = ndarray.sum(result, 0).view(self.dtype)

        # Return result
        return result
    
    def display(self, time = slice(None), reactant_normalized = False):
        """
        display returns a string that represents the net
        reaction.  The net reaction can be optionally sliced
        by time and divided by the largest reactant stoichiometry
        """
        value = self[time].sum()        
        positive = value >= 0
        
        if positive:
            return "%f gained by %s" % (value, self.name)
        else:
            return "%f lost by %s" % (value, self.name)

    def array(self):
        """
        view IPR array
        """
        return self.view(type = ndarray, dtype = 'f')
        
