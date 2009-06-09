__all__ = ['Species', 'Stoic', 'Reaction', 'Process', 'IPR', 'Mechanism', \
           'mechanisms', 'getmech', 'get_mech', 'get_pure_mech', \
           'get_prepared_mech', 'analyses']

if __name__ != '__main__':
    from SpeciesGroup import Species
    from ReactionGroup import Reaction, Stoic
    from ProcessGroup import Process
    from IPRArray import IPR
    from Mechanism import Mechanism
    import mechanisms
    import getmech
    from getmech import get_pure_mech, get_prepared_mech
    get_mech = get_pure_mech
else:
    from optparse import OptionParser
    from types import MethodType
    from warnings import warn

    parser = OptionParser()
    parser.set_usage("Usage: %prog [-c cb05_camx|cbiv_camx] [mrgfile]")
    parser.add_option("-i", "--interactive", dest="interactive", \
                        action="store_true", default=False, \
                        help="open an interactive environment", \
                        metavar="interactive environment")
    
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

        from analyses.net_balance import net_balance
        dict(net_balance=net_balance)[options.analysis](args, options)
        
    from perm import mechanisms, \
                            netcdf, \
                            getmech
    
    get_prepared_mech = getmech.get_prepared_mech
    get_pure_mech = getmech.get_pure_mech
    try:
        mech = get_prepared_mech(options.mechanism)
    except:
        mech = get_pure_mech(options.mechanism)
        
    if len(args) > 0:
        NetCDFFile = netcdf.NetCDFFile
        mrg_file = NetCDFFile(args[0],'rs')
        mech.set_mrg(mrg_file)

    def load_environ(locals_dict):
        if not locals_dict.has_key('mech'):
            locals_dict['mech'] = mech
        locals_dict.update(mech.species_dict)
        locals_dict.update(mech.reaction_dict)
        try:
            locals_dict.update(mech.nreaction_dict)
        except:
            pass
    
        locals_dict.update([(k,getattr(mech,k)) for k in dir(mech) if '__' not in k and isinstance(getattr(mech,k),MethodType) and k not in ('set_mrg', 'set_irr', 'set_ipr')])            
        
    if options.interactive:
        from perm.Shell import PERMConsole
        console = PERMConsole()
        load_environ(console.locals)
        console.interact()
    else:
        load_environ(globals())
        for script in args[1:]:
            execfile(script)

