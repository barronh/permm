__all__ = ['Species', 'Stoic', 'Reaction', 'Process', 'PhyTables', \
           'PtbTable', 'NetTables', 'SumTable', 'IPR', 'Mechanism', \
           'mechanisms', 'IRRTable', 'VOCTable', 'getmech', \
           'get_mech', 'get_pure_mech', 'get_prepared_mech']

if __name__ != '__main__':
    from SpeciesGroup import Species
    from ReactionGroup import Reaction, Stoic
    from ProcessGroup import Process
    from PhyTableMaker import PhyTable
    from IPRArray import IPR
    from Mechanism import Mechanism
    from SumTableMaker import SumTable
    from NetTableMaker import NetTables
    from PhyTableMaker import PhyTables, VOCTable
    from PtbMaker import PtbTable
    from IRRTable import IRRTable
    import mechanisms
    from getmech import get_pure_mech, get_prepared_mech
    get_mech = get_pure_mech
else:
    from net_yaml import *
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-c cb05_camx|cbiv_camx] [mrgfile]")
    parser.add_option("-i", "--interactive", dest="interactive", \
                        action="store_true", default=False, \
                        help="open an interactive environment", \
                        metavar="interactive environment")
    
    parser.add_option("-c", "--mechanism", dest="mechanism", \
                      default="cb05_camx", help="Chemical mechanisms", \
                      metavar="MECHANISM")

    parser.add_option("-o", "--output", dest="output", \
                      default=None, help="Base output path", \
                      metavar="ouput base")
    
    (options, args) = parser.parse_args()
    if len(args)<1:
        if options.interactive:
            mrg_data_path = None
        else:
            parser.error(msg="Requires a pyPA mrg file as an argument for analysis output.  You can enter an interactive environment (-i), but will not have access to \"netted\" reactions")
    else:
        mrg_data_path=args[0]

    from net_yaml.mechanisms.cb05_camx.mechprep import cb05_camx_prep as mech_prep
    from net_yaml.mechanisms.getmech import get_prepared_mech, get_pure_mech
    
    #mech_yaml_path = '/Users/barronh/Development/net_yaml/mechanisms/cb05_camx/cb05_camx.yaml'
    mech = get_prepared_mech(options.mechanism)
    if mrg_data_path is not None:
        from net_yaml.netcdf import NetCDFFile
        mrg_file = NetCDFFile(mrg_data_path,'rs')
        mech.set_mrg(mrg_file)
        
    if options.interactive:
        from warnings import warn
        globals().update(mech.species_dict)
        globals().update(mech.reaction_dict)
        try:
            globals().update(mech.nreaction_dict)
        except:
            pass

        while False:
            try:
                cmd = raw_input('>> ')
                exec(cmd)
            except:
                warn('Use quit() to exit')
    else:
        if options.output is None:
            net_data_path = '.'.join(mrg_data_path.split('.')[:-1])+'.net.'
        else:
            net_data_path = options.output+'.'
        print >> file(net_data_path+'sum','wb'), SumTable(mech)
        print >> file(net_data_path+'net','wb'), NetTables(mech)
        print >> file(net_data_path+'phy','wb'), PhyTables(mech)
        print >> file(net_data_path+'voc','wb'), VOCTable(mech)
        print >> file(net_data_path+'ptb','wb'), PtbTable(mech)
        mech = get_pure_mech(options.mechanism)
        if mrg_data_path is not None:
            from pynetcdf import NetCDFFile
            mrg_file = NetCDFFile(mrg_data_path,'rs')
            mech.set_mrg(mrg_file)
        
        print >> file(net_data_path+'irr','wb'), IRRTable(mech)