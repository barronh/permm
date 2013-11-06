def MakeMultiDiGraph(mech, traceable = None, slice = slice(None)):
    import networkx as nx
    from numpy import sum
    G = nx.MultiDiGraph()
    for spc in traceable:
        for rxn in mech.get_irrs(reactants = spc):
            rxn_slice = rxn[slice]
            rcts = [mech(rct) for rct in rxn.reactants()]
            rct_frac = rxn_slice[spc]/sum([rxn_slice[rct] for rct in rcts], axis = 0)
            for prd_name in rxn.products():
                prd = mech(prd_name)
                if prd in traceable:
                    G.add_edge(spc, prd, (rxn_slice[prd] * rct_frac).sum())
                    
    return G

def MakeDiGraph(mech, traceable = None, slice = slice(None)):
    import networkx as nx
    G = MakeMultiDiGraph(mech = mech, traceable = traceable, slice = slice)
    NewG = nx.DiGraph()
    for edge in G.edges_iter():
        u,v = edge
        vals = G[u][v]
        NewG.add_edge(u, v, sum(vals))
        
    return NewG
    
def MakeUniDiGraph(mech, traceable = None, slice = slice(None)):
    import networkx as nx
    G = MakeDiGraph(mech = mech, traceable = traceable, slice = slice)
    NewG = nx.DiGraph()
    g_edges = G.edges()
    for edge in g_edges:
        u,v = edge
        data = G.get_edge_data(u, v)
        if NewG.has_edge(v, u):
            new_data = NewG.get_edge_data(v, u) - data
            if data > 0:
                NewG[v][u] = new_data
            else:
                NewG.delete_edge(v, u)
                NewG.add_edge(u, v, abs(new_data))
        else:
            NewG.add_edge(u, v, data)
    
    return NewG

def MakeCarbonTrace(mech, slice = slice(None), traceable = None, makes_larger = []):
    from numpy.ma import masked_invalid
    traceable = [spc for spcn, spc in mech.species_dict.iteritems() if any([spcv.get(spc.name, {}).get('C', 0) > 0 for spcv in  spc.spc_dict.values()])]
    print traceable
    #if len(traceable) == 0: traceable = mech('[eval(n) for n in HC.names()]')
    import networkx as nx
    from numpy import sum
    G = nx.MultiDiGraph()
    for spc in traceable:
        spc_C = spc.atoms('C')[spc]
        for rxn in mech.get_irrs(reactants = spc):
            rxn_slice = rxn[slice].sum()
            rcts = [mech(rct) for rct in rxn.reactants()]
            for prd_name in rxn.products():
                prd = mech(prd_name)
                try:
                    prd_C = prd.atoms('C')[prd]
                except KeyError:
                    # Product has no carbon
                    continue

                if prd_C <= spc_C or spc in makes_larger:
                    if prd in traceable:
                        rct_frac = rxn_slice[spc]/sum([rxn_slice[rct] for rct in rcts if rct in traceable and rct.atoms('C')[rct] >= prd_C], axis = 0)
                        spc_portion = masked_invalid(rxn_slice[prd] * rct_frac).sum()
                        print spc, spc_C, prd, prd_C, rct_frac
                        G.add_edge(spc, prd, weight = spc_portion)
                    
    return G
