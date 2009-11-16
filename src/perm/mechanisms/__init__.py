__all__ = ['cb05_camx', 'cb05_cmaq', 'geos_chem', 'saprc99_cmaq', 'saprc07_cmaq', 'racm2_cmaq', 'leedsmcm_kpp', 'mz4_kpp']
from ..getmech import get_pure_mech

for mech in __all__:
    exec('def %s():\n    return get_pure_mech("%s")' % (mech, mech), globals(), locals())
