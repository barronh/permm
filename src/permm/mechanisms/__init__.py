__all__ = ['cb05_camx', 'cb05_cmaq', 'geos_chem', 'saprc99_cmaq', 'saprc07_cmaq', 'racm2_cmaq', 'mz4_kpp', 'small_strato']
from ..getmech import get_pure_mech
from yaml import load
from os.path import dirname, \
                    abspath, \
                    join
atoms_path = abspath(join(dirname(__file__), 'atoms.yaml'))

atoms = load(file(atoms_path).read())
for mech in __all__:
    exec('def %s():\n    return get_pure_mech("%s")' % (mech, mech), globals(), locals())
