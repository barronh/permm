__all__ = ['Mechanism', 'Species', 'Reaction', 'atoms', 'mechanism_dict']

from .core.Species import Species
from .core.Reaction import Reaction
from .core.Mechanism import Mechanism
from .mechanisms import mechanism_dict, atoms
from . import mechanisms
from . import getmech
from . import Shell
from . import GUI
from .getmech import get_pure_mech, get_prepared_mech
get_mech = get_pure_mech

if __name__ == '__main__':
    from permm.main import parse_and_run
    parse_and_run()