__all__ = ['Mechanism', 'Species', 'Reaction', 'atoms', 'mechanism_dict']

if __name__ != '__main__':
    from core.Species import Species
    from core.Reaction import Reaction
    from core.Mechanism import Mechanism
    from mechanisms import mechanism_dict, atoms
    import mechanisms
    import getmech
    import Shell
    import GUI
    from getmech import get_pure_mech, get_prepared_mech
    get_mech = get_pure_mech
else:
    from permm.main import parse_and_run
    parse_and_run()