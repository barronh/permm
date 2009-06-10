__all__ = ['Species', 'Stoic', 'Reaction', 'Process', 'IPR', 'Mechanism', \
           'mechanisms', 'getmech', 'get_mech', 'get_pure_mech', \
           'get_prepared_mech', 'analyses', 'graphing', 'GUI', 'Shell']

if __name__ != '__main__':
    from SpeciesGroup import Species
    from ReactionGroup import Reaction, Stoic
    from ProcessGroup import Process
    from IPRArray import IPR
    from Mechanism import Mechanism
    import mechanisms
    import getmech
    import Shell
    import GUI
    from getmech import get_pure_mech, get_prepared_mech
    get_mech = get_pure_mech
else:
    from optparse import OptionParser
    from warnings import warn
    from perm.mechanisms import __all__ as all_mechs
    from perm.analyses import __all__ as all_analyses
    all_mechs = '|'.join(all_mechs)
    all_analyses = '|'.join(all_analyses)
    parser = OptionParser()
    parser.set_usage("Usage: python -m perm [-c %s] [-i|-a %s] [mrgfile]" % (all_mechs, all_analyses))
    parser.add_option("-i", "--interactive", dest="interactive", \
                        action="store_true", default=False, \
                        help="open an interactive environment", \
                        metavar="interactive environment")
    
    parser.add_option("-g", "--gui", dest="graphical", \
                        action="store_true", default=False, \
                        help="open a graphical user interactive environment", \
                        metavar="graphical user interactive environment")
    
    parser.add_option("-c", "--mechanism", dest="mechanism", \
                      default="cb05_camx", help="Chemical mechanisms", \
                      metavar="MECHANISM")

    parser.add_option("-a", "--analysis", dest="analysis", \
                      default=None, help="Stock analysis to perform", \
                      metavar="ANALYSIS")

    parser.add_option("-o", "--output", dest="output", \
                      default=None, help="Base output path", \
                      metavar="ouput base")
    
    (options, args) = parser.parse_args()

    if options.analysis is not None:
        if len(args)<1:
            parser.error(msg="Requires a pyPA mrg file as an argument for analysis output.  You can enter an interactive environment (-i), but will not have access to \"netted\" reactions")

        from perm.analyses.net_balance import net_balance
        dict(net_balance=net_balance)[options.analysis](args, options)
        
    from perm import mechanisms, \
                            netcdf, \
                            getmech
    from perm.Shell import load_environ
    
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech
    mech = get_pure_mech(options.mechanism)
        
    if len(args) > 0:
        NetCDFFile = netcdf.NetCDFFile
        mrg_file = NetCDFFile(args[0],'rs')
        mech.set_mrg(mrg_file)

        
    load_environ(mech, globals())
    for script in args[1:]:
        execfile(script)

    if options.graphical:
        from perm.GUI import StartGUI
        StartGUI(mech)
    elif options.interactive:
        from perm.Shell import PERMConsole
        console = PERMConsole()
        load_environ(mech,console.locals)
        console.interact()
    elif len(args) == 0:
        parser.print_usage()
