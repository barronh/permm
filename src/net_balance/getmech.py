from net_yaml import Mechanism

def get_pure_mech(mechanism):
    from os.path import dirname, \
                        abspath, \
                        join
    mech_path = abspath(join(dirname(__file__),mechanism,'%s.yaml' % (mechanism,)))
    return Mechanism(mech_path)
    
def get_prepared_mech(mechanism):
    from net_yaml.mechanisms.cb05_camx.mechprep import cb05_camx_prep as cb05_camx
    from net_yaml.mechanisms.cb05_cmaq.mechprep import cb05_cmaq_prep as cb05_cmaq
    from net_yaml.mechanisms.geos_chem.mechprep import geos_chem_prep as geos_chem
    mech_object = get_pure_mech(mechanism)
    mech_prep = eval(mechanism)

    return mech_prep(mech_object)