__all__ = ['mechprep', 'mech']

import mechprep
from ...getmech import get_pure_mech
from os import path
mech = get_pure_mech(path.basename(path.dirname(__file__)))
