__all__ = ['mz4_chem_prep']

def z4_chem_prep(mech_object):
    from permm.Mechanism import Mechanism
    from permm.SpeciesGroup import Species

    globals().update(mech_object.species_dict)
    
    return mech_object
