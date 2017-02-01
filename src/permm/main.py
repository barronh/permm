def parse_and_run():
    import os
    from argparse import ArgumentParser
    from warnings import warn
    from permm import mechanism_dict, Mechanism
    from permm.analyses import __all__ as all_analyses
    all_mechs = '|'.join(list(mechanism_dict.keys()))
    all_analyses = '|'.join(all_analyses)
    from PseudoNetCDF.pncparse import getparser, pncparse
    parser = ArgumentParser(description = "permm (Python Environment for Reaction Mechanism Mathematics)")
    parser.add_argument("--gui", dest="graphical", \
                        action="store_true", default=False, \
                        help="open a graphical user interactive environment")
    
    parser.add_argument("--mechanism", dest="mechanism", \
                      default="cb05_camx", help="Chemical mechanisms (e.g., a custom path|%s)" % all_mechs)

    parser.add_argument("--analysis", dest="analysis", \
                      default=None, help="Stock analysis to perform (i.e., %s)" % all_analyses)
    parser.add_argument("--scripts", dest = "scripts", action = 'append', default = [], help = "Provide scripts to run")
    parser.add_argument("-i", "--interactive", dest = "interactive", action = 'store_true', default = False, help = "Run interactively")
    parser.add_argument("--input-role", dest = "inrole", action = 'append', default = [], help = "Is this a 'merge' file or a 'concentration' file or any other word indicating type?")
    
    parser.add_argument("pseudonetcdf", nargs = '*')
    options = parser.parse_args()
    if len(options.pseudonetcdf) > 0:
        pncparser = getparser(has_ofile = False, interactive = True)
        ifiles, pncoptions = pncparse(has_ofile = False, interactive = True, parser = pncparser, args = args.pseudonetcdf)
    else:
        ifiles = []
        pncoptions = {}

    from permm import netcdf
    from permm.Shell import load_environ
    try:
        mech = mechanism_dict[options.mechanism]
    except KeyError:
        mech = Mechanism(options.mechanism)
        
    start_script = 0
    options.inrole += ['merge'] * (len(ifiles) - len(options.inrole))
    for r, ifile in zip(options.inrole, ifiles):
        if r == 'merge':
            mech.set_mrg(ifile)
        else:
            vardict = dict([(k, v) for k, v in list(ifile.variables.items()) if k in mech.species_dict])
            mech.set_process(r, vardict)
            
    from permm.Shell import PERMConsole
    console = PERMConsole()
    load_environ(mech,console.locals)

    for script in options.scripts:
        if os.path.isfile(script):
            exec(compile(open(script).read(), script, 'exec'), globals(), console.locals)
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
            console.runsource("history = matrix(mech, [C2O3], [HC+Radical-OH-HO2-O1D], [])")
            console.runsource("history.run()")
        else:
            raise "Unkown analysis"

    if options.interactive:
        load_environ(mech,console.locals)
        console.interact()

        
if __name__ == '__main__':
    parse_and_run()