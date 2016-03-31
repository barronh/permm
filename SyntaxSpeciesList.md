# Species Syntax #

Each species contains a label and an atom-based definition.  The label or species name can be any YAML recognizable string (e.g. `NO2` or `'NO'`).  The atom-based definition can combine multiples of any atoms using the pattern below.  If you would prefer not to define a species, simply type IGNORE.

Examples:
```
    'NO': 1*N + 1*O
    ALK4: 4*C
    O3: 3*O
    CL: 1*Cl
    NR: IGNORE
```

Valid Patterns:
```
    label: number*atom + number*atom + ...
```

or

```
    label: IGNORE
```

