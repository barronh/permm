__all__ = ['Species', 'Stoic', 'Reaction', 'Process', 'PhyTables', \
           'PtbTable', 'NetTables', 'SumTable', 'IPR', 'Mechanism', \
           'mechanisms', 'IRRTable', 'VOCTable', 'getmech', \
           'get_mech', 'get_pure_mech', 'get_prepared_mech']

if __name__ != '__main__':
    from SpeciesGroup import Species
    from ReactionGroup import Reaction, Stoic
    from ProcessGroup import Process
    from IPRArray import IPR
    from Mechanism import Mechanism
    from PhyTableMaker import PhyTable
    from SumTableMaker import SumTable
    from NetTableMaker import NetTables
    from PhyTableMaker import PhyTables, VOCTable
    from PtbMaker import PtbTable
    from IRRTable import IRRTable
    import mechanisms
    import getmech
    from getmech import get_pure_mech, get_prepared_mech
    get_mech = get_pure_mech
else:
    from optparse import OptionParser
    from types import MethodType
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

    from net_balance import mechanisms, \
                            netcdf, \
                            getmech
    
    from net_balance.PhyTableMaker import PhyTable
    from net_balance.SumTableMaker import SumTable
    from net_balance.NetTableMaker import NetTables
    from net_balance.PhyTableMaker import PhyTables
    from net_balance.PhyTableMaker import VOCTable
    from net_balance.PtbMaker import PtbTable
    from net_balance.IRRTable import IRRTable

    mech_prep = mechanisms.cb05_camx.mechprep.cb05_camx_prep
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech
    try:
        mech = get_prepared_mech(options.mechanism)
    except:
        mech = get_pure_mech(options.mechanism)
        
    if mrg_data_path is not None:
        NetCDFFile = netcdf.NetCDFFile
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
        
        globals().update([(k,getattr(mech,k)) for k in dir(mech) if '__' not in k and isinstance(getattr(mech,k),MethodType) and k not in ('set_mrg', 'set_irr', 'set_ipr')])

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