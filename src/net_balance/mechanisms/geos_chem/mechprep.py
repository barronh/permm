__all__ = ['geos_chem_prep']

def geos_chem_prep(mech_object):
    from net_yaml.Mechanism import Mechanism
    from net_yaml.SpeciesGroup import Species

    globals().update(mech_object.species_dict)
    
    return mech_object
