import code
import readline
import atexit
import os
from types import MethodType

class PERMConsole(code.InteractiveConsole):
    def __init__(self, locals = None, filename = '<console>',
                histfile = os.path.expanduser("~/.permm-history")):
        code.InteractiveConsole.__init__(self)
        self.init_history(histfile)
    
    def init_history(self, histfile):
        readline.parse_and_bind("tab: complete")
        if hasattr(readline, "read_history_file"):
            try:
                readline.read_history_file(histfile)
                readline.set_history_length(1000)
            except IOError:
                pass
            atexit.register(self.save_history, histfile)
            
    def save_history(self, histfile):
        readline.write_history_file(histfile)

def load_environ(mech, locals_dict):
    if not locals_dict.has_key('mech'):
        locals_dict['mech'] = mech
    locals_dict.update(mech.species_dict)
    locals_dict.update(mech.reaction_dict)
    try:
        locals_dict.update(mech.irr_dict)
        locals_dict.update(mech.nreaction_dict)
    except:
        pass
    try:
        locals_dict.update(mech.process_dict)
    except:
        pass

    locals_dict.update([(k,getattr(mech,k)) for k in dir(mech) if '__' not in k and isinstance(getattr(mech,k),MethodType) and k not in ('set_mrg', 'set_irr', 'set_ipr')])
    