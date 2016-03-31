# Tutorial 1: Use the gui to create a plot of the top producers of acetaldehyde #
With the graphical user interface, we are going to query the mechanism and create a publication quality graphic.

Step 1: Open a terminal emulator

Open any old terminal emulator.  Navigate to your test data.  If you need help, ask a friend.

Step 2: Load data and start application

```
$ python -m permm -g test.mrg.nc
```
You should have gotten a little window.  

Step 3: Select the reactants

For any query, we select the reactants using the "Select Reactants" list.  This time we don't care what the reactants are.  In fact, we don't need to know.  Leave it blank.

Step 4: Set the query and/or

For any query, we select if want reactions with the selected reactants AND products or the reactants OR products.  In this case, we want the default AND.  Leaving reactants blank was like saying "any reactant", so if we said "any reactant OR ALD2 products"... we'd get any reaction.

Step 3: Select the products

For any query, we select the products using the "Select Products" list.  This tutorial is for acetaldehyde, so I'm choosing ALD2.  You can choose another compound if you'd like.

Step 5: Plot the Reactions

We want to make a plot, so we choose "Plot Reactions" and hit "go"

Step 6: Review and save

Look at your graph and, if you like it, click the save button.  Save it as a PDF to get the best publication quality.

Step 7:

Close the gui.