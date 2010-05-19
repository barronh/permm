from Tkinter import Checkbutton, Frame, Label, Scrollbar, Listbox, Button, IntVar, Tk, VERTICAL, EXTENDED, END, N, S
import os
from types import MethodType

class TkApp:
    def __init__(self, master, mech):
        self.mech = mech
        frame = Frame(master)
        frame.grid()
        rct_label = Label(frame, text = 'Select Reactants')
        rct_label.grid(column = 1, row = 1)
        and_or = Label(frame, text = 'AND/OR')
        and_or.grid(column = 3, row = 1)
        prod_label = Label(frame, text = 'Select Products')
        prod_label.grid(column = 4, row = 1)
        what_to_do = Label(frame, text = 'Execute')
        what_to_do.grid(column = 6, row = 1)
        reactants_scrollbar = Scrollbar(frame, orient = VERTICAL)
        self.reactants = Listbox(frame, selectmode = EXTENDED, exportselection = 0, yscrollcommand = reactants_scrollbar.set)
        self.reactants.grid(column = 1, row = 2)
        reactants_scrollbar.config(command = self.reactants.yview)
        reactants_scrollbar.grid(column = 2, row = 2, sticky = N + S)

        self.logical_and = IntVar()
        self.logical_and.set(1)
        c = Checkbutton(frame, text = "AND", variable = self.logical_and)
        c.grid(column = 3, row = 2)
        
        products_scrollbar = Scrollbar(frame, orient = VERTICAL)
        self.products = Listbox(frame, selectmode = EXTENDED, exportselection = 0)
        self.products.grid(column = 4, row = 2)
        products_scrollbar.config(command = self.products.yview)
        products_scrollbar.grid(column = 5, row = 2, sticky = N + S)
        
        self.method_list = Listbox(frame, selectmode = EXTENDED, exportselection = 0)
        self.method_list.grid(column = 6, row = 2)
        #self.methods = [k for k in dir(self.mech) if k[:1] != '_' and isinstance(getattr(self.mech, k), MethodType)]
        self.methods = ['plot_rxns', 'find_rxns', 'print_rxns', 'print_irrs', 'print_net_rxn', 'plot_proc']
        method_labels= ['Plot Reactions', 'Show Rxn Ids', 'Print Rxns', 'Print IRRs', 'Print Net Rxn', 'Process Plot']
        for method in method_labels:
            self.method_list.insert(END, method)
            
        species_keys = self.mech.species_dict.keys()
        species_keys.sort()
        self.species_objects = [self.mech.species_dict[spc] for spc in species_keys]

        for spc in species_keys:
            self.reactants.insert(END, spc)
            self.products.insert(END, spc)
        self.execute_button = Button(frame, text = "go", command = self.execute)
        self.execute_button.grid(column = 6, row = 4)

    def _get_species(self, list):
        items = list.curselection()
        try: items = map(int, items)
        except: pass
        items = map(lambda i, d = self.species_objects: d[i], items)
        return items
        
    def get_reactants(self):
        return self._get_species(self.reactants)
        
    def get_products(self):
        return self._get_species(self.products)
        
    def get_methods(self):
        items = self.method_list.curselection()
        try: items = map(int, items)
        except: pass
        items = map(lambda i, d = self.methods: d[i], items)
        return items
        
    def execute(self):
        os.system('clear')
        reactants = self.get_reactants()
        products = self.get_products()
        logical_and = bool(self.logical_and.get())
        methods = self.get_methods()
        kwds = dict(reactants = reactants, products = products, logical_and = logical_and)
        for method in methods:
            if method == 'plot_proc':
                thiskwds = {}
                thiskwds['species'] = (reactants+products)[0]
            else:
                thiskwds = kwds.copy()
            getattr(self.mech, method)(**thiskwds)

def StartTk(mech):
    root = Tk()
    
    app = TkApp(root, mech)
    
    root.mainloop()
