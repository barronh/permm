from numpy import array, zeros
def cycles(mech, species_group, indirect = False):
    mech.globalize(globals())
    species_group_names = species_group.names()
    species_group_cycles = {}
    for spc_name in species_group.names():
        spc = eval(spc_name)
        spc_names = set(spc.names())
        this_spc_cycles = species_group_cycles[spc] = {}
        init = make_net_rxn(-species_group, spc)[spc]
        this_spc_cycles['source'] = init
        prop_total = zeros(init.shape, dtype = 'f')
        for prop in get_irrs(species_group, spc):
            made = prop.get_spc(spc)
            prop_total += made
            rad_rcts = [eval(rct) for rct in set(prop.reactants()).intersection(species_group_names)]
            rad_rct_stoics = array([prop.get_spc(rad, 'r') for rad in rad_rcts])
            rad_rct_stoic = rad_rct_stoics.sum(0)
            rad_frac = rad_rct_stoics / rad_rct_stoic[None,:]
            for radi, rad in enumerate(rad_rcts):
                if (made != 0).any():
                    this_spc_cycles.setdefault('propfrom', {})[rad] = made / rad_frac[radi]

        sink_total = zeros(init.shape, dtype = 'f')
        template = zeros(init.shape, dtype = 'f')
        for sink in get_irrs(spc, -species_group):
            made = sink.get_spc(spc)
            sink_total += made
            nrad_rcts = [eval(rct) for rct in set(sink.reactants()).difference(spc_names)]
            if nrad_rcts == []:
                nrad_rcts = [spc]
            nrad_rct_stoics = array([sink.get_spc(rad, 'r') for rad in nrad_rcts])
            nrad_rct_stoic = nrad_rct_stoics.sum(0)
            nrad_frac = nrad_rct_stoics / nrad_rct_stoic[None,:]
            for nradi, nrad in enumerate(nrad_rcts):
                if (made != 0).any():
                    this_spc_cycles.setdefault('sinkby', {}).setdefault(nrad, template.copy())[:] = made / nrad_frac[nradi]
        
        this_spc_cycles['reacted'] = make_net_rxn(spc)[spc]
        this_spc_cycles['prop'] = prop_total.copy()
        this_spc_cycles['sink'] = sink_total.copy()

        exec('propeff = 1 - sink / reacted', globals(), this_spc_cycles)
        exec('total = source + prop', globals(), this_spc_cycles)
        
    if indirect:
        for spc in this_spc_cycles.keys():
            this_spc_cycles = species_group_cycles[spc]
            for prop_source, prop in this_spc_cycles['propby'].iteritems():
                propeff = species_group_cycles[prop_source]['propeff']
                
    return species_group_cycles
    