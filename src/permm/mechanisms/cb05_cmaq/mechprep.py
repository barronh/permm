__all__ = ['cb05_cmaq_prep']

def cb05_cmaq_prep(mech_object):
    from ...Mechanism import Mechanism
    from ...SpeciesGroup import Species
    from pynetcdf import NetCDFFile

    globals().update(mech_object.species_dict)
    
    ### Decorate mechanism with pseudo-Species
    # RO2 place holder
    mech_object.add_spc_to_reactions(mech_object.species_dict['RO2'])    
    
    # Secondary aldehydes
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['FORM'], name = 'sFORM'))
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['ALD'], name = 'sALD'))
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['ALD_'], name = 'sALD_'))
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['ALD2'], name = 'sALD2'))
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['ALDX'], name = 'sALDX'))
    mech_object.add_prd_to_reactions(Species(mech_object.species_dict['MGLY'], name = 'sMGLY'))
    mech_object.species_dict['sFORM'] = Species(names = ['sFORM'], stoic = [1.], name = 'sFORM')
    mech_object.species_dict['sALD'] = Species(names = ['sALD'], stoic = [1.], name = 'sALD')
    mech_object.species_dict['sALD_'] = Species(names = ['sALD_'], stoic = [1.], name = 'sALD_')
    mech_object.species_dict['sALD2'] = Species(names = ['sALD2'], stoic = [1.], name = 'sALD2')
    mech_object.species_dict['sALDX'] = Species(names = ['sALDX'], stoic = [1.], name = 'sALDX')
    mech_object.species_dict['sMGLY'] = Species(names = ['sMGLY'], stoic = [1.], name = 'sMGLY')

    # Add xHO2 when it is only a fraction of HO2 is secondary
    xHO2 = Species(name = 'xHO2', names = ['XO2'], stoic = [1.], exclude = False)
    cond = lambda x: x.has_prd(XO2) and x.has_prd(HO2) and x.get(XO2,0) <= x.get(HO2, 0)
    mech_object.add_spc_to_reactions_with_condition(xHO2,cond)
    
    # Add xHO2 when all HO2 is secondary
    xHO2 = Species(name = 'xHO2', names = ['HO2'], stoic = [1.], exclude = False)
    cond = lambda x: x.has_prd(XO2) and x.has_prd(HO2) and x.get(XO2,0) > x.get(HO2, 0)
    mech_object.add_spc_to_reactions_with_condition(xHO2,cond)
    
    # Add xHO2 when all HO2 is and HO2 production is explicit
    cond = lambda x: x.has_rct(XO2_N_TO2) and x.has_prd(HO2)
    mech_object.add_spc_to_reactions_with_condition(xHO2,cond)
    
    # Add primar HO2
    pHO2 = Species(name = 'pHO2', names = ['HO2', 'xHO2'], stoic = [1., -1])
    mech_object.add_prd_to_reactions(pHO2)
    
    # Add consumed HO2
    lHO2 = Species(name = 'lHO2', names = ['HO2'], stoic = [1.])
    mech_object.add_rct_to_reactions(lHO2)

    # Over-write default instance of lHO2
    lHO2 = Species(name = 'lHO2', names = ['lHO2'], stoic = [1.])
    mech_object.species_dict['lHO2'] = lHO2

    # Over-write default instance of pHO2
    pHO2 = Species(name = 'pHO2', names = ['pHO2'], stoic = [1.])
    mech_object.species_dict['pHO2'] = pHO2

    # Over-write default instance of xHO2
    xHO2 = Species(name = 'xHO2', names = ['xHO2'], stoic = [1.])
    mech_object.species_dict['xHO2'] = xHO2

    return mech_object

