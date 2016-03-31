# Tutorial 3: Use a script to create a plot of the top producers of acetaldehyde #
With a script, we are going to query the mechanism and create a publication quality graphic.

Step 1: Open your favorite text editor

Step 2: Type the following code
```
plot_rxns(products = ALD2, path = 'ALD2.pdf')
```

Step 3: Save the script
Save the script as ald2.py in the same place as the test data you downloaded

Step 4: Open a terminal emulator
Open any old terminal emulator.  Navigate to the test data.  If you need help, ask a friend.

Step 5: Load data and run the script
```
$ python -m permm -g test.mrg.nc ald2.py
```
