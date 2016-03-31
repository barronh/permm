# Plotting the NOx/VOC Sensitivity #

## Introduction ##
A system is NOx or VOC sensitive when adding or removing NOx or VOC changes the O3 concentration.  Sandy Sillman developed an indicator of sensitivity based on the reactions of radicals.  If NOx is in abundance (i.e. VOC-sensitive), radicals react with NOx and form HNO3.  If NOx is rare, radicals react with other radicals to form peroxides.  By examining the ratio of peroxide to HNO3 production, we can classify an air parcel as NOx or VOC-sensitive (aka limited).

## Assumption ##
This script will assume you are using the reaction mechanism Carbon Bond '05, but it is easily adaptable to any mechanism.

## Details ##
With a script, we are going to perform 4 operations.  First, we are going to define peroxides.  Second, we are going to create a net reaction for HNO3 production.  Third, we are going to create a net reaction for peroxides.  Fourth, we are going to create a ratio.  Fifth, we are going to plot the result.

## Steps ##
Step 1: Open your favorite text editor

Step 2: Add the definition of peroxides
```
PEROXIDES = ROOH + H2O2
```

Step 3: Add the following net reaction definitions
```
nrxn_prod_peroxides = make_net_rxn(products = PEROXIDES)
nrxn_prod_hno3 = make_net_rxn(products = HNO3)
```

Step 4: Create the ratios
```
prooh = nrxn_prod_peroxides[PEROXIDES]
phno3 = nrxn_prod_hno3[HNO3]
sillman_ratio = prooh/phno3
```

Step 5: Plot the ratio
```
plot(sillman_ratio, path = 'sillman_ratio.pdf')
```

Step 6:
Save the script as sillman_ratio.py in the same place as the test data you downloaded

Step 7: Open a terminal emulator
Open any old terminal emulator.  Navigate to the test data.  If you need help, ask a friend.

Step 8: Load data and run the script
```
$ python -m permm -g test.mrg.nc sillman_ratio.py
```
