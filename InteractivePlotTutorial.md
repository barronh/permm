# Tutorial 2: Use the terminal to create a plot of the top producers of acetaldehyde #
With the terminal, we are going to query the mechanism and create a publication quality graphic.

Step 1: Open a terminal emulator
Open any old terminal emulator.  If you need help, ask a friend.

Step 2: Load data and start application
```
$ python -m permm -i -c cb05_camx test.mrg.nc
```
You should have gotten a little window.  

Step 3: Review the syntax for plot_rxns
```
>>> plot_rxns(products = ALD2)
```

Step 4: Review and save
Look at your graph and, if you like it, click the save button.  Save it as a PDF to get the best publication quality.

Step 5: Close the graph and exit the terminal
Use the x to close the graph and then type exit() to close the terminal