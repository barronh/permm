__all__=['timeseries', 'barplot']
if __name__ != '__main__':
    import timeseries
    import barplot
else:
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-t] [-m <model>]  <yamlfile>")
    parser.add_option("-t", "--template", dest="template",action="store_true",default=False,
                        help="Output template on standard out (configurable with -m and -c", metavar="Template")
    
    parser.add_option("-m", "--model", dest="model",default="new",
                        help="Model can either be camx, cmaq or wrfchem (for use with -t)", metavar="MODEL")
    
    (options, args) = parser.parse_args()
    if options.template:
        from os.path import join,abspath,dirname
        import sys
        cmaq_template_path = join(abspath(dirname(__file__)),'defaults','template.yaml')       
        print >> sys.stdout, file(cmaq_template_path,'r').read()
        parser.exit()
    if len(args)<1:
        parser.error(msg="Requires a yaml file as an argument.  For a template use the -t option.  The template will be output to the stdout.")
        parser.exit()
    else:
        import yaml
        from permm.graphing.timeseries import phy_plots, chem_plots, rxn_plots
        conf = yaml.load(file(args[0]))
        phy_plots(conf)
        chem_plots(conf)
        rxn_plots(conf)
