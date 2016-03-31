# Mechanism Syntax #
## Overview ##
PERMM Mechanisms are a collection of chemical species definitions, reaction definitions and species group definitions.  These definitions are saved in a YAML file that shares syntax with KPP, the LaTeX package mhchem, and attempts to follow traditional chemistry syntax.  Each chemical species is a collection of atoms.  Each species group is a collection of chemical species. Each reaction combines species with reaction roles (reactant/product) and stoichiometry.  

 
## Example Mechanism Definition Yaml ##
```
---
species_list:
    CL: Cl
    O1D: O
    'NO': N + O
    ALK4: 4*C
    ...
reaction_list:
    IRR_1:  O3 + NO ->[k] 1.000*NO2
    IRR_2:  NO2 ->[j] 1.000*NO + 1.000*O
    IRR_3:  O2 + O ->[k] O3
    ...

species_group_list:
- NOx = NO + NO2
- NOz = HNO3 + RNO3
- PANS = PAN + HPAN + GPAN
- NOy = NOx + 2*N2O5 + NO3 + HNO4 + NOz + PANS
...
```

more information on each list are available on subpages:
 - [species_list](SyntaxSpeciesList.md)
 - [reaction_list](SyntaxReactionList.md)
 - [species_group_list](SyntaxSpeciesGroupList.md)
