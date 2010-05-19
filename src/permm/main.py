def parse_and_run():
    import os
    from optparse import OptionParser
    from warnings import warn
    from permm.mechanisms import __all__ as all_mechs
    from permm.analyses import __all__ as all_analyses
    all_mechs = '|'.join(all_mechs)
    all_analyses = '|'.join(all_analyses)
    parser = OptionParser()
    parser.set_usage("Usage: python -m permm [-c %s] [-g|-i|-a %s] [mrgfile]" % (all_mechs, all_analyses))
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

    from permm import mechanisms, \
                            netcdf, \
                            getmech
    from permm.Shell import load_environ
    
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech

    mech = get_pure_mech(options.mechanism)
    start_script = 0
    if len(args) > 0:
        try:
            try:
                NetCDFFile = netcdf.NetCDFFile
                mrg_file = NetCDFFile(args[0],'r')
            except:
                from PseudoNetCDF.camxfiles.Memmaps import irr
                NetCDFFile = irr
                mrg_file = NetCDFFile(args[0])

            mech.set_mrg(mrg_file)
           
            
        except (IOError, RuntimeError), (e):
            warn("\n-------------------------------------------------------\nFirst argument was not a data file.\nAttempting to continue with first argument as a script.\n-------------------------------------------------------")
        except ValueError, (e):
            warn("\n-------------------------------------------\nAre you sure this is a %s data file?\n-------------------------------------------" % options.mechanism)
            raise e
        else:
            start_script = 1

    from permm.Shell import PERMConsole
    console = PERMConsole()
    load_environ(mech,console.locals)

    for script in args[start_script:]:
        if os.path.isfile(script):
            execfile(script, globals(), console.locals)
        else:
            exec(script, globals(), console.locals)


    if options.graphical:
        console.runsource("from permm.GUI import StartGUI")
        console.runsource("StartGUI(mech)")

    if options.analysis is not None:
        if len(args)<1:
            parser.error(msg="Requires a pyPA mrg file as an argument for analysis output.  You can enter an interactive environment (-i), but will not have access to \"netted\" reactions")

        if options.analysis == "net_balance":
            console.runsource("from permm.analyses.net_balance import net_balance")
            console.runsource("net_balance('%s', '%s', '%s')" % (options.mechanism, args[0], options.output))
        elif options.analysis == "history":
            console.runsource("from permm.analyses.history import matrix")
            console.runsource("history = matrix(mech, [C2O3], [HC], [])")
            console.runsource("history.run()")
        else:
            raise "Unkown analysis"

    if options.interactive:
        load_environ(mech,console.locals)
        console.interact()

        
if __name__ == '__main__':
    parse_and_run()