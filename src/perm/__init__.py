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
    import os
    from optparse import OptionParser
    from warnings import warn
    from perm.mechanisms import __all__ as all_mechs
    from perm.analyses import __all__ as all_analyses
    all_mechs = '|'.join(all_mechs)
    all_analyses = '|'.join(all_analyses)
    parser = OptionParser()
    parser.set_usage("Usage: python -m perm [-c %s] [-g|-i|-a %s] [mrgfile]" % (all_mechs, all_analyses))
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

    from perm import mechanisms, \
                            netcdf, \
                            getmech
    from perm.Shell import load_environ
    
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech

    mech = get_pure_mech(options.mechanism)
    
    if len(args) > 0:
        try:
            NetCDFFile = netcdf.NetCDFFile
            mrg_file = NetCDFFile(args[0],'rs')
            mech.set_mrg(mrg_file)
        except:
            start_script = 0
        else:
            start_script = 1

    from perm.Shell import PERMConsole
    console = PERMConsole()
    load_environ(mech,console.locals)

    def runsource(source, filename = '<input>'):
        for line in source.splitlines():
            console.runsource(line, filename, 'single')
            
    for script in args[start_script:]:
        if os.path.isfile(script):
            source = file(script).read()
            fname = script
        else:
            source = script
            fname = '<input>'

        runsource(source, filename = fname)

    if options.interactive:
        load_environ(mech,console.locals)
        console.interact()

    if options.graphical:
        console.runsource("from perm.GUI import StartGUI")
        console.runsource("StartGUI(mech)")

    if options.analysis is not None:
        if len(args)<1:
            parser.error(msg="Requires a pyPA mrg file as an argument for analysis output.  You can enter an interactive environment (-i), but will not have access to \"netted\" reactions")

        if options.analysis == "net_balance":
            console.runsource("from perm.analyses.net_balance import net_balance")
            console.runsource("net_balance('%s', '%s', '%s')" % (options.mechanism, args[0], options.output))
        
