from pdb import set_trace
from numpy import zeros, exp, newaxis

class matrix():
    def __init__(self, mech, node, traceable, intermediates, bottomup = True):
        self.__mech = mech
        self.__traceable = set(traceable)
        self.__intermediates = set(intermediates)
        self.traced = set([node])
        self.origins = {}
        self.origin_loss = {}
        self.producers = {}
        self.concentrations = {}
        self.production = {}
        self.losses = {}
        self.__ntimes = mech.mrg.variables['TFLAG'][:,0,0].shape[0]
        self.intermediates = set(intermediates)
        if bottomup:
            self.add_predecessors(node)
        else:
            self.add_successors(node)
        self.update_origin_contributions()
    
    def add_predecessors(self, node):
        mech = self.__mech
        producers = mech.find_rxns(products = mech(node))
        losers = mech.find_rxns(reactants = mech(node))

        ztemp = zeros(self.__ntimes, dtype = 'f')
        try:
            HT = mech('H_Trans[%s].array()' % (node,))
            VT = mech('V_Trans[%s].array()' % (node,))
            PE = mech('Emissions[%s].array()' % (node,))
            LDep = mech('Deposit[%s].array()' % (node,))

            PHT = HT.copy(); PHT[HT < 0] = 0
            LHT = HT.copy(); LHT[HT > 0] = 0
            PVT = VT.copy(); PVT[VT < 0] = 0
            LVT = VT.copy(); LVT[VT > 0] = 0
        except KeyError:
            PHT = LHT = PVT = LVT = PE = LDep = ztemp.copy()
        
        try:
            self.concentrations.setdefault(node,{})['Initial'] = mech('Initial[%s].array()' % (node,))
            self.concentrations[node]['Final'] = mech('Final[%s].array()' % (node,))
            average = self.concentrations[node]['Average'] = mech('0.5 * (Final[%s].array()+Initial[%s].array())' % (node,node))
        except:
            self.concentrations.setdefault(node,{})['Initial'] = ztemp.copy()
            self.concentrations[node]['Final'] = ztemp.copy()
            
        
        self.producers[node] = dict(HT = PHT,
                                     VT = PVT,
                                     E = PE,
                                     C = ztemp.copy())
        self.losses[node] = LHT + LVT + LDep
        
        total_loss = self.losses[node]
        for rxn in losers:
            total_loss += mech(rxn)[node]
        

        for rxn in producers:
            reactants = mech.reaction_dict[rxn].reactants()
            produces = mech(rxn)[node]
            self.producers[node]['C'] += produces
            tracers = self.__traceable.intersection(reactants)
            ntracers = len(tracers)
            if ntracers > 1:
                produces = produces/ntracers
            
            for tracer in tracers:
                node_producers = self.producers.setdefault(node,{})
                if node_producers.has_key(tracer):
                    node_producers[tracer] += produces
                else:
                    node_producers[tracer] = produces
                
                if tracer not in self.traced:
                    self.traced.add(tracer)
                    self.add_predecessors(tracer)

        self.production[node] = PHT + PVT + PE + self.producers[node]['C']

        if not self.concentrations[node].has_key('Average'):
            average = self.concentrations[node]['Average'] = 0.5 * self.production[node]

        total_loss /= average

    def update_origin_contributions(self):
        mech = self.__mech
        for node, origin_productions in self.producers.iteritems():
            origin_contributions = self.origins[node] = {}
            loss_contributions = self.origin_loss[node] = {}
            origin_contributions['Initial'] = zeros(self.__ntimes + 1, dtype = 'f')
            origin_contributions['Initial'][0] = self.concentrations[node]['Initial'][0]
            
            hourly_loss = self.losses[node]
            for origin, hourly_production in origin_productions.iteritems():
                hourly_origin_contribution = origin_contributions[origin] = zeros(self.__ntimes + 1, dtype = 'f')
                hourly_loss_contribution = loss_contributions[origin] = zeros(self.__ntimes + 1, dtype = 'f')
                for old_origin_contribution, new_origin_contribution, loss_contribution, loss, prod in zip(hourly_origin_contribution[:-1][:,newaxis], hourly_origin_contribution[1:][:,newaxis], hourly_loss_contribution[1:][:,newaxis], hourly_loss, hourly_production):
                    if origin == 'C' and False:
                        set_trace()
                    new_origin_contribution.itemset(prod/-loss * (1 - exp(loss)) +  old_origin_contribution * exp(loss))
                    loss_contribution.itemset(old_origin_contribution + prod - new_origin_contribution)
                    
