# Reaction Syntax #

Reactions are combinations of species, stoichiometry, reaction role, and reaction type.  Each reaction is identified by a reaction label, which can be any YAML recognizable string that starts with a letter and contains only letters and numbers.  The reaction definition starts with a reactants list containing one or more reactants with optional stoichiometry.  The reactants list is followed by an arrow (i.e. `->`) with the reaction type enclosed in braces (i.e. `[k]` or `[j]`) where j specifies photolysis and k specifies thermal.  The reaction arrow and type are followed by an optional products list that has the same syntax as the reactants list.

Reaction definition regular expression:
```(?P<reactants>.*)=(?P<rxn_type>[kj])[>]\s*(?P<products>.*)```

Reactant/product regular expression:
```(\s?(?P<sign>[+-])?\s?)?((?P<stoic>\d(\.(\d{1,3}(E\d{2})?)?)?)\*)?(?P<name>[A-Z]\w{0,3})\s?```

Example definition:
```
    IRR_1: OH + OLE =k> 0.8*FORM + 0.33*ALD2 + 0.62*ALDX + 0.8*XO2 + 0.95*HO2 - 0.7 PAR
```
