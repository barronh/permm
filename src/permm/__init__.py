__all__ = ['Species', 'Stoic', 'Reaction', 'Process', 'IPR', 'Mechanism', \
           'mechanisms', 'getmech', 'get_mech', 'get_pure_mech', \
           'get_prepared_mech', 'analyses', 'graphing', 'GUI', 'Shell', 'guis']

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
    from permm.main import parse_and_run
    parse_and_run()