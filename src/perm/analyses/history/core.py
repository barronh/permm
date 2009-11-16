from pdb import set_trace
from numpy import zeros, exp, newaxis, rollaxis
from perm.SpeciesGroup import Species

class matrix(object):
    def __init__(self, mech, nodes, traceable, intermediates, bottomup = True):
        """
            mechanism - perm.Mechanism.Mechanism object
            nodes - list of species to trace
            intermediates - species that should be immediately attributed to predecessors
            bottomup - trace from node to root? no, then root to node
        """
        # Internalize arguments for later use
        self.__mech = mech
        self.__traceable = set(traceable)
        self.__intermediates = set(intermediates)
        self.__bottomup = bottomup
        
        # make sure that nodes is a list of nodes
        if isinstance(nodes, Species):
            nodes = [nodes]

        # Initialize dictionaries for storing
        # source data
        self.origins = {}
        self.origin_loss = {}

        # Producers is a dictionary (keys = species; values = dictionaries of producers) 
        # that contains integrated production arrays
        #
        # Example: {'CO': {'Transport': array([8,5,6,...]), ...}}
        self.producers = {}
        
        # Concentrations is a dictionary (keys = species; values dictionary 
        # initial, average, and final concentration arrays.  Similar to producers
        self.concentrations = {}
        
        # Production is a dictionary of gross production each species
        self.production = {}
        
        # Losses is a dictionary of gross consumption of each species
        self.losses = {}

        # Preset traced nodes as those provided by the user
        self.traced = set(nodes)
        
        
        # Storing data about input shapes
        self.__old_shape = mech('INIT').shape
        self.__shape = [i for i in self.__old_shape]
        self.__ntimes = mech.mrg.variables['TFLAG'].shape[-3]
        for i, s in enumerate(self.__shape):
            if s == self.__ntimes:
                self.__shape[i] = self.__ntimes + 1
                self.__time_dim = i
        self.__shape = tuple(self.__shape)
        
        self.run()
        
    def run(self):
        if self.__bottomup:
            for node in [n for n in self.traced]:
                self.add_predecessors(node)
        else:
            for node in [n for n in self.traced]:
                self.add_successors(node)

        self.update_origin_contributions()
    
    def add_predecessors(self, node):
        """
            Find species that lead to the production of node
        """
        
        mech = self.__mech
        
        # Find all reactions that produce node
        producers = mech.find_rxns(products = node)
        
        # Find all reactions that consume
        consumers = mech.find_rxns(reactants = node)

        # Create a temporary array
        ztemp = zeros(self.__old_shape, dtype = 'd')

        try:
            HT = mech('H_Trans')[node].array()
            VT = mech('V_Trans')[node].array()
            M = mech('Motion')[node].array()
            PE = mech('Emissions')[node].array()
            LDep = mech('Deposit')[node].array()

            PHT = HT.copy(); PHT[HT < 0] = 0
            LHT = HT.copy(); LHT[HT > 0] = 0
            PVT = VT.copy(); PVT[VT < 0] = 0
            LVT = VT.copy(); LVT[VT > 0] = 0
            PM = M.copy(); PM[VT < 0] = 0
            LM = M.copy(); LM[VT > 0] = 0
        except (KeyError, NameError), (e):
            PM = LM = PHT = LHT = PVT = LVT = PE = LDep = ztemp.copy()
        
        # Add initial, final, and (if available average) concentration
        try:
            self.concentrations.setdefault(node.name,{})['Initial'] = mech('Initial')[node].array()
            self.concentrations[node.name]['Final'] = mech('Final')[node].array()
            average = self.concentrations[node.name]['Average'] = 0.5 * (mech('Final')[node].array() + mech('Initial')[node].array())
        except:
            self.concentrations.setdefault(node.name,{})['Initial'] = ztemp.copy()
            self.concentrations[node.name]['Final'] = ztemp.copy()
            
        
        # Populate producers dictionary for node
        self.producers[node.name] = dict(HT = PHT,
                                     VT = PVT,
                                     E = PE,
                                     M = PM,
                                     C = ztemp.copy())
        
        # Populate losses dictionary for node with physical
        # losses
        total_loss = self.losses[node.name] = LHT + LVT + LDep + LM
        
        # Add chemical losses
        for rxn in consumers:
            total_loss += mech(rxn)[node]
        
        # For each producing reaction, trace traceable reactants
        for rxn in producers:
            # Get reactants of this reaction
            reactants = mech.reaction_dict[rxn].reactants()
            
            # Get total net production of node by this reaction
            produces = mech(rxn)[node]
            
            # Increment chemistry process
            self.producers[node.name]['C'] += produces
            
            # Create a list of traceable species in reaction reactants
            tracers = reduce(list.__add__, [list(set(tr.names()).intersection(reactants)) for tr in self.__traceable])
            ntracers = len(tracers)
            
            # The reactions production of node, must be divided among n tracers
            if ntracers > 1:
                produces = produces/ntracers

            # Create output for clarity
            print mech.reaction_dict[rxn]
            print " -",

            # For each traceable species in reactants
            for tracer in tracers:
                print "%s," % tracer,
                node_producers = self.producers.setdefault(node.name,{})
                if node_producers.has_key(tracer):
                    node_producers[tracer] += produces
                else:
                    node_producers[tracer] = produces
            print ""
            
            for tracer in tracers:
                if mech(tracer) not in self.traced:
                    print "Tracing", tracer, "from", mech.reaction_dict[rxn]
                    self.traced.add(mech(tracer))
                    self.add_predecessors(mech(tracer))
        self.production[node.name] = PHT + PVT + PE + self.producers[node.name]['C']

        if not self.concentrations[node.name].has_key('Average'):
            average = self.concentrations[node.name]['Average'] = 0.5 * self.production[node.name]

        total_loss /= average

    def update_origin_contributions(self):
        """
        """
        mech = self.__mech
        
        for node, origin_productions in self.producers.iteritems():
            # Capture the origins of node by source
            origin_sources = self.origins[node] = {}
            # and the loss of origins by source
            origin_losses = self.origin_loss[node] = {}

            # Load sources with initial values
            origin_sources['Initial'] = zeros(self.__shape, dtype = 'd')
            origin_sources['Initial'][..., :-1] = self.concentrations[node]['Initial'][..., :]
            
            # Losses are not origin specific
            hourly_loss = self.losses[node]
            hourly_loss = rollaxis(hourly_loss, self.__time_dim)
            
            for origin, hourly_production in origin_productions.iteritems():
                # Create an array to house production fro this origin
                # store it in the self.origins dictionary
                origin_sources[origin] = hourly_origin_sources = zeros(self.__shape, dtype = 'd')
                origin_losses[origin] = hourly_origin_losses = zeros(self.__shape, dtype = 'd')
                
                # Creating pointer arrays with the hour dimension first
                # for ease of iteration
                hourly_origin_sources = rollaxis(hourly_origin_sources,self.__time_dim)
                hourly_origin_losses = rollaxis(hourly_origin_losses, self.__time_dim)
                hourly_production = rollaxis(hourly_production, self.__time_dim)
                
                for thishr in range(1, self.__ntimes+1):
                    ## Create pointers to loss and production rates
                    loss = hourly_loss[thishr-1]
                    prod = hourly_production[thishr-1]

                    ## Create pointers to the origin arrays
                    # the last hour's origin
                    old_origin = hourly_origin_sources[thishr-1]
                    # this hour's origin
                    new_origin = hourly_origin_sources[thishr]
                    # this hour's origin loss
                    origin_loss = hourly_origin_losses[thishr]

                    # Debugging entry point                    
                    if origin == 'C' and False:
                        set_trace()
                    
                    # Capture production and loss terms based on
                    # an ordinary differential equation approximation
                    # See Tonnesen 1995 Dissertation
                    this_prod = prod/-loss * (1 - exp(loss)) +  old_origin * exp(loss)
                    this_loss = old_origin + prod - new_origin
                    
                    # Update origin arrays; testing for scalar
                    try:
                        new_origin.itemset(this_contribution)
                        origin_loss.itemset(this_loss)
                    except:
                        new_origin[...] = this_prod
                        origin_loss[...] = this_loss

    def traceback(self, node, start = -1):
        for spc in node.names():
            pass
        raise
            