from pdb import set_trace
from numpy import zeros, exp, newaxis, rollaxis
from permm.core.Species import Species

class matrix(object):
    def __init__(self, mech, nodes, traceable, intermediates, bottomup = True):
        """
            mechanism - permm.Mechanism.Mechanism object
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
        tmp_init = mech('INIT')
        self.__old_shape = tmp_init[tmp_init.keys()[0]].shape
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
    
    def add_predecessors(self, node, debug = False):
        """
        Recursively find reactions that consume and/or produce chemical 
        species (node)_n and attribute production to traceable reactant 
        species (node)_{n-1}.
        """
        
        mech = self.__mech
        
        # Find all reactions that produce node
        producers = mech.find_rxns(products = node)
        
        # Find all reactions that consume
        consumers = mech.find_rxns(reactants = node)

        # Create a temporary array
        ztemp = zeros(self.__old_shape, dtype = 'd')

        try:
            Initial = mech('Initial')[node].array()
            H_Trans = mech('H_Trans')[node].array()
            V_Trans = mech('V_Trans')[node].array()
            Motion = mech('Motion')[node].array()
            Emissions = mech('Emissions')[node].array()
            Deposition = mech('Deposit')[node].array()
            
            # Initial concentrations are production terms only for the first hour
            rollaxis(Initial, self.__time_dim)[1:] = 0
            ProdH_Trans = H_Trans.copy(); ProdH_Trans[H_Trans < 0] = 0
            LossH_Trans = H_Trans.copy(); LossH_Trans[H_Trans > 0] = 0
            ProdV_Trans = V_Trans.copy(); LossV_Trans[V_Trans < 0] = 0
            LossV_Trans = V_Trans.copy(); LossV_Trans[V_Trans > 0] = 0
            ProdMotion = Motion.copy(); ProdMotion[Motion < 0] = 0
            LossMotion = Motion.copy(); LossMotion[Motion > 0] = 0
            
            # Populate producers dictionary for node
            self.producers[node.name] = dict(H_Trans = ProdH_Trans,
                                             V_Trans = ProdV_Trans,
                                             Emissions = Emissions,
                                             Motion = ProdMotion,
                                             Initial = Initial,
                                             Chemistry = ztemp.copy())
            # Populate losses dictionary for node with physical
            # losses
            total_loss = self.losses[node.name] = LossH_Trans + LossV_Trans + Deposition + LossMotion
        except (KeyError, NameError), (e):
            # Populate producers dictionary for node with zeros
            self.producers[node.name] = dict(H_Trans = ztemp.copy(),
                                             V_Trans = ztemp.copy(),
                                             Emissions = ztemp.copy(),
                                             Motion = ztemp.copy(),
                                             Chemistry = ztemp.copy(),
                                             Initial = ztemp.copy())
            # Populate losses dictionary for node with physical
            # losses
            total_loss = self.losses[node.name] = ztemp.copy()
        
        # Add initial, final, and (if available average) concentration
        try:
            self.concentrations.setdefault(node.name,{})['Initial'] = mech('Initial')[node].array()
            self.concentrations[node.name]['Final'] = mech('Final')[node].array()
            average = self.concentrations[node.name]['Average'] = 0.5 * (mech('Final')[node].array() + mech('Initial')[node].array())
        except:
            self.concentrations.setdefault(node.name,{})['Initial'] = ztemp.copy()
            self.concentrations[node.name]['Final'] = ztemp.copy()
            
                
        
        # Add chemical losses
        for rxn in consumers:
            total_loss += mech(rxn)[node]
        
        # For each producing reaction, trace traceable reactants
        for rxn in producers:
            # Get reactants of this reaction
            reactants = mech.reaction_dict[rxn].reactants()
            
            # Get total net production of node by this reaction
            produces = mech(rxn)[node].copy()
            
            # Increment chemistry process
            self.producers[node.name]['Chemistry'] += produces
            
            # Create a list of traceable species in reaction reactants
            tracers = reduce(list.__add__, [list(set(tr.names()).intersection(reactants)) for tr in self.__traceable])
            ntracers = len(tracers)
            
            # The reactions production of node, must be divided among n tracers
            if ntracers > 1:
                produces = produces/ntracers

            # Create output for clarity
            if debug:
                print mech.reaction_dict[rxn]
                print " -",

            # For each traceable species in reactants
            for tracer in tracers:
                if debug:
                    print "%s," % tracer,
                    
                # Attribute tracer production of node to tracer
                node_producers = self.producers.setdefault(node.name,{})
                if node_producers.has_key(tracer):
                    node_producers[tracer] += produces
                else:
                    node_producers[tracer] = produces
            
            if debug:
                print ""
            
            # For each traceable reactant
            for tracer in tracers:
                # Check that it has not already been traced
                if mech(tracer) not in self.traced:
                    if debug:
                        print "Tracing", tracer, "from", mech.reaction_dict[rxn]
                    
                    # Store tracer as having been traced
                    self.traced.add(mech(tracer))

                    # Trace species
                    self.add_predecessors(mech(tracer))

        self.production[node.name] = sum([prod for prod in self.producers[node.name].values()])
        
        # Total 1st order loss rate is approximated by total 
        # loss divided by average concentration
        if not self.concentrations[node.name].has_key('Average'):
            average = self.concentrations[node.name]['Average'] = eval('.5 * (Initial + Final)', None, self.concentrations[node.name])
        else:
            average = self.concentrations[node.name]['Average']

        total_loss /= average

    def update_origin_contributions(self):
        """
        
        """
        mech = self.__mech
        
        for node, producers in self.producers.iteritems():
            # Capture the origins of node by source
            origin_sources = self.origins[node] = {}
            # and the loss of origins by source
            origin_losses = self.origin_loss[node] = {}

            # Load sources with initial values
            origin_sources['Initial'] = zeros(self.__shape, dtype = 'd')
            origin_sources['Initial'][..., :-1] = self.concentrations[node]['Initial'][..., :]
            
            # Losses are not origin specific
            loss_rate = self.losses[node]
            loss_rate = rollaxis(loss_rate, self.__time_dim)
            
            for origin, production in producers.iteritems():
                # Create an array to house production fro this origin
                # store it in the self.origins dictionary
                origin_sources[origin] = sources = zeros(self.__shape, dtype = 'd')
                origin_losses[origin] = losses = zeros(self.__shape, dtype = 'd')
                
                # Creating pointer arrays with the hour dimension first
                # for ease of iteration
                sources = rollaxis(sources,self.__time_dim)
                losses = rollaxis(losses, self.__time_dim)
                production = rollaxis(production, self.__time_dim)
                
                for thishr in range(1, self.__ntimes+1):
                    ## Create pointers to loss and production rates
                    loss = loss_rate[thishr-1]
                    prod = production[thishr-1]

                    ## Create pointers to the origin arrays
                    # the last hour's origin
                    old_origin = sources[thishr-1]

                    # Debugging entry point                    
                    if origin == 'C' and False:
                        set_trace()
                    
                    # Capture production and loss terms based on
                    # an ordinary differential equation approximation
                    # See Tonnesen 1995 Dissertation
                    sources[thishr] = new_origin = prod/-loss * (1 - exp(loss)) +  old_origin * exp(loss)
                    losses[thishr] = old_origin + prod - new_origin
    
    def reattribute_origins(self):
        """
        Timestep contributions from each chemical species
        should be attributed to the fractional origin of that species
        for each timestep.
        """
        fractional_origins = {}
        
        for ti in range(self.__ntimes):
            for created_spc, origins in self.origins.keys():
                while any([_k in self.__mech.species_dict for _k in origins.keys()]):
                    for origin_key in origins.keys():
                        if origin_key in self.__mech.species_dict:
                            origin_origins = self.origins.get(origin_key, {})
                                    
                    