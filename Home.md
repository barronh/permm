# Welcome to PERMM #

PERMM is a Python-based Environment for Reaction Mechanisms/Mathematics.  Simply put, PERMM helps you detangle the complex relationships in reaction mechanisms.  Reaction mechanisms, as shown in our logo, often have complex and recursive relationships that when combined with 4-dimensional rate data... can get frizzy. PERMM starts with a reaction mechanism, applies [n-dimensional rate data](ProcessAnalysisData.md), and then provides an interface for chemical analysis.

PERMM lets you:
  * define chemical species families
  * query reactions by species (or family)
  * net reactions against each other
  * perform math over n-dimensions
  * and so much more

Below is a series of examples commands to illustrate the simplicity and power of PERMM.  The first example will return all reaction names that create an "Odd Oxygen" from a radical.  The second example will print them all.  The third will create a net reaction for a time-series, sum it across time and display the integrated net reaction.  Each function uses the same arguments, and this continuity helps make PERMM easy to learn.

```
>>> find_rxns(Radical, Ox)
['IRR_10', 'IRR_103', 'IRR_107', 'IRR_108', 'IRR_132', 'IRR_26', 'IRR_28', 'IRR_29', 'IRR_30', 'IRR_31', 'IRR_33', 
  'IRR_41', 'IRR_47', 'IRR_48', 'IRR_54', 'IRR_61', 'IRR_68', 'IRR_81', 'IRR_88', 'IRR_89', 'IRR_92']
>>> print_rxns(Radical, Ox)
IRR_10 1.0*O1D =k> 1.0*O
IRR_103 1.0*NO + 1.0*CXO3 =k> 1.0*ALD2 + 1.0*HO2 + 1.0*NO2 + 1.0*XO2
IRR_107 1.0*PANX + 1.0*OH =k> 1.0*ALD2 + 1.0*NO2
IRR_108 1.0*HO2 + 1.0*CXO3 =k> 0.2*AACD + 0.2*O3 + 0.8*PACD
IRR_132 1.0*TO2 + 1.0*NO =k> 0.1*NTR + 0.9*HO2 + 0.9*NO2 + 0.9*OPEN
etc...
>>> make_net_rxn(Radical, Ox).sum()
0.9*NTR + 1.0*TO2 + 1.0*PANX + 1.0*O1D + 1.0*NO3 + 1.0*HONO + 1.0*HCO3 + 2.0*CXO3 + 3.0*C2O3 + 7.0*NO + 8.0*OH =k> 
  -0.66*PAR + 0.33*ALDX + 0.4*AACD + 0.4*O3 + 0.9*HO2 + 0.9*OPEN + 1.0*FACD + 1.0*PAN + 1.33*FORM + 1.6*PACD + 2.0*HNO3 + 
   2.0*O + 2.33*ALD2 + 7.9*NO2
```

PERMM can be used with scripts, interactively, or as a graphical user interface.  To learn about your options, use the help option.

```
$ python -m permm --help
```

PERMM provides many commonly used atmospheric reaction mechanisms as implemented by their host models, but you can also define your own.  First, check to see whether it already exists.  We already have CBIV for CAMx, CB05 for CMAQ and CAMx, SAPRC99 for CMAQ, SAPRC07 for CMAQ, and GEOS-Chem v08-01-01.  Run PERMM with the help option to see a list of supported mechanisms.  If we don't have your chemical mechanism, you can [define your own](MechanismSyntax.md).

First, you need to [get PERMM installed](InstallPERMM.md).  To get a feel for the basic functionality, use the gui.  To learn more about using PERMM with rate data, see our [tutorial](Tutorial.md).

Enjoy!<br />
[The Development Team](Development.md)

p.s. Our logo is a networkx diagram of an atmospheric oxidation reaction mechanism.