import code
import readline
import atexit
import os

class PERMConsole(code.InteractiveConsole):
    def __init__(self, locals = None, filename = '<console>',
                histfile = os.path.expanduser("~/.perm-history")):
        code.InteractiveConsole.__init__(self)
        self.init_history(histfile)
    
    def init_history(self, histfile):
        readline.parse_and_bind("tab: complete")
        if hasattr(readline, "read_history_file"):
            try:
                readline.read_history_file(histfile)
            except IOError:
                pass
            atexit.register(self.save_history, histfile)
            
    def save_history(self, histfile):
        readline.write_history_file(histfile)
