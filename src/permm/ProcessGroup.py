from utils import AttrDict
__all__ = ['Process']

class Process(AttrDict):
    """
    Process is a group of processes where the simplest case is
    a 1 process group where names[0] == name
    
    The easiest way to create a process is with keywords:
        Process(
                name = '<process name>',
                names = ['component name1', ..., 'component name2']
               )
        
    Processes can be added together to create larger processes
        >>> A_W
        {'exclude': False, 'names': ['A_W'], 'name': 'A_W'}
        >>> A_E
        {'exclude': False, 'names': ['A_E'], 'name': 'A_E'}
        >>> H_Trans = A_W + A_E + A_S + A_N + D_S + D_W + D_N + D_E
        >>> H_Trans.name = 'H_Trans'
        >>> H_Trans
        {
            'name': 'H_Trans',
            'names': ['A_W', 'A_E', 'A_S', 'A_N', 'D_S', 'D_W', 'D_N', 'D_E'],
            'exclude': False
        }
    """
    def __init__(self, *args, **kwds):
        """
        Processes require a name (string) and a names (list)
        """
        AttrDict.__init__(self, *args, **kwds)
        
        has_prc_name = self.has_key('name')
        has_prc_names = self.has_key('names')
        
        self.setdefault('exclude', False)
        
        if not (has_prc_name and has_prc_names):
            raise KeyError, "ProcessGroupDict must be initialized with name and names"

    def __add__(self,y):
        """
        Process1 + Process2 is a group process containing the 
        components of each
        """
        if isinstance(y,Process):
            if y.exclude:
                new_names = list(set(self.names).difference(y.names))
            else:
                new_names = list(set(self.names).union(y.names))
            new_name = self.name + '+' + y.name
        else:
            raise TypeError, "Only processes can be added together"
            
        return Process(name =  new_name, names = new_names)