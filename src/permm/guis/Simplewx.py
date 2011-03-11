import sys
import wx
from permm.Shell import PERMConsole, load_environ

class WXAgg(wx.Panel):
    def __init__(self, parent, mech):
        wx.Panel.__init__(self, parent)
        self.mech = mech
        self.species = mech.species_dict.keys()
        self.species.sort()

        # A button
        self.gobutton =wx.Button(self, label="Go", pos=(200, 325))
        self.Bind(wx.EVT_BUTTON, self.OnGo, self.gobutton)

        # the combobox Control
        self.methods = ['plot_rxns', 'find_rxns', 'print_rxns', 'print_irrs', 'print_net_rxn', 'plot_proc', 'make_net_rxn']
        self.method_labels= ['Plot Reactions', 'Show Rxn Ids', 'Print Rxns', 'Print IRRs', 'Print Net Rxn', 'Process Plot', 'Print Net Reaction']
        self.method_select = wx.ComboBox(self, size=(-1, -1), choices=self.method_labels, style=wx.CB_READONLY)
        self.method_select.SetValue(self.method_labels[0])

        self.reactants_select = wx.ListBox(self, size = (-1, 200), choices = self.species, style = wx.LB_MULTIPLE)
        self.products_select = wx.ListBox(self, size = (-1, 200), choices = self.species, style = wx.LB_MULTIPLE)
 
        # Checkbox
        self.logical_and = wx.CheckBox(self, label="And?", pos=(20,180))
        self.logical_and.SetValue(True)
        
        sizer = wx.FlexGridSizer(rows=4, cols=3, hgap=10, vgap=5)
        sizer.AddGrowableCol(3)
        sizer.AddSpacer(1) # R1, C1
        sizer.Add(wx.StaticText(self, label="PERMM")) # R1, C2
        sizer.AddSpacer(1) # R1, C3
        sizer.Add(wx.StaticText(self, label="Reactants")) # R2, C1
        sizer.AddSpacer(1) # R2, C2
        sizer.Add(wx.StaticText(self, label="Products")) # R2, C3
        sizer.Add(self.reactants_select) # R3, C1
        sizer.Add(self.logical_and) # R3, C2
        sizer.Add(self.products_select) # R3, C3
        sizer.AddSpacer(self.method_select) # R4, C1
        sizer.AddSpacer(1) # R4, C2
        sizer.Add(self.gobutton) # R4, C3

        sizer1 = wx.FlexGridSizer(rows=1, cols=2, hgap=10, vgap=5)
        # A multiline TextCtrl - This is here to show how the events work in this program, don't pay too much attention to it
        self.logger = wx.TextCtrl(self, size=(400,300), style=wx.TE_MULTILINE | wx.TE_READONLY)
        sizer1.Add(sizer)
        sizer1.Add(self.logger)
        
        self.SetSizer(sizer1)
        self.Fit()

        



    def GetReactants(self):
        ids = self.reactants_select.GetSelections()
        spcs = [self.species[id] for id in ids]
        return spcs

    def GetProducts(self):
        ids = self.products_select.GetSelections()
        spcs = [self.species[id] for id in ids]
        return spcs
    
    def GetMethod(self):
        return self.methods[self.method_labels.index(self.method_select.GetValue())]

    def OnGo(self,event):
        method = self.GetMethod()
        logical_and = self.logical_and.IsChecked()
        reactants = ', '.join(self.GetReactants())
        products = ', '.join(self.GetProducts())
        eval_str = '%s(reactants = [%s], products = [%s], logical_and = %s)' % (method, reactants, products, logical_and)
        self.logger.AppendText(eval_str + ':\n')
        return_text = self.mech(eval_str)
        if return_text is not None:
            print return_text

def StartWx(mech):
    app = wx.App(True)
    frame = wx.Frame(None)
    panel = WXAgg(frame, mech)
    frame.Fit()
    frame.Show()
    app.MainLoop()

if __name__ == '__main__':
    from permm.mechanisms import cb05_cmaq
    mech = cb05_cmaq()
    StartWx(mech)