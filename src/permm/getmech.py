from core.Mechanism import Mechanism

def get_pure_mech(mechanism):
    from os.path import dirname, \
                        abspath, \
                        join, \
                        exists
    if exists(mechanism) and '.yaml' in mechanism:
        mech_path = mechanism
    else:
        mech_path = abspath(join(dirname(__file__),'mechanisms', '%s.yaml' % (mechanism,)))
    if not exists(mech_path):
        raise ImportError, "Mechanism you supplied is not a known mechanism and is not a file path to a definition"
        
    return Mechanism(mech_path)
    
def get_prepared_mech(mechanism):
    from mechanisms.cb05_camx.mechprep import cb05_camx_prep as cb05_camx
    from mechanisms.cb05_cmaq.mechprep import cb05_cmaq_prep as cb05_cmaq
    from mechanisms.geos_chem.mechprep import geos_chem_prep as geos_chem
    mech_object = get_pure_mech(mechanism)
    mech_prep = eval(mechanism)

    return mech_prep(mech_object)