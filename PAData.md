# Process Analysis Data #

Process Analysis (PA) data is a tally of the mass of chemical species that go through some process.  The processes available are defined by the operator splitting of model, but typically include advection, diffusion, deposition, emissions, aerosol processing, and chemistry.  Reactions are an expanded view of chemistry.  Process Analysis data keeps track of the total mass that goes through each process and each reaction.

```
NOx_FINAL = NOx_INIT + SUM(NOx_PROCESSES)

NOx_PROCESSES = NOx_ADVECTION + NOx_DIFFUSION + NOx_DEPOSITION + ...

NOx_ADVECTION = NOx_WEST_ADV + NOx_EAST_ADV + ... + NOx_TOP_ADV + NOx_BOTTOM_ADV

NOx_CHEM = NOx_PROD + NOx_LOSS
```