import re
import sys
import yaml

from Mechanism import Mechanism

reaction_pattern = re.compile("^(?P<r_num>\d{1,3})\s(?P<reactants>[A-Z]\w{0,3}(\s\+\s[A-Z]\w{0,3})*)\s\n?(?P<products>(\d(\.(\d{1,3}(E\d{2})?)?)?\s)?[A-Z]\w{0,3}(\s\n?(\+|\-)\n?\s\n?(\d(\.(\d{1,3}(E\d{2})?)?)?\s\n?)?[A-Z]\w{0,3})*)?\s\n?(?P<reaction_type>Photolysis|(\d\.\d{3}E[+-]\d{2}))",re.M)

species_pattern = re.compile("(\s?\n?(?P<sign>[+-])?\s?\n?)?(?P<stoic>\d(\.(\d{1,3}(E\d{2})?)?)?\s\n?)?(?P<name>[A-Z]\w{0,3})\s?",re.M)
reactions_text = file("../../mechanisms/cb05_camx/cb05_camx_pdf.txt","r").read()

reaction_list = [match.groupdict() for match in reaction_pattern.finditer(reactions_text)]

tab = "    "
exclude_species = ("O2", "M", "H2O", "None")
species_list = []
reaction_yaml = "reaction_list:\n"
species_yaml = "species_list:\n"

def consolidate_dups(spc_group):
    spc_group_names = set([s['name'] for s in spc_group])
    if len(spc_group_names) < len(spc_group):
        new_spc_group=[dict(name=name) for name in spc_group_names]
        for new_spc in new_spc_group:
            for old_spc in [old_spc for old_spc in spc_group if old_spc['name'] == new_spc['name']]:
                if new_spc.get('sign',old_spc['sign']) != old_spc['sign']:
                    raise ValueError, 'Math not implemented for different signs'
                else:
                    new_spc['stoic'] = str(float(new_spc.get('stoic',0))+float(old_spc['stoic']))
                    new_spc['sign'] = old_spc['sign']
        return new_spc_group
    else:
        return spc_group

def remove_nones(spc_group):
    spc_group_names = [s['name'] for s in spc_group]
    while None in spc_group_names:
        spc_group.pop(spc_group_names.index(None))
        spc_group_names.pop(spc_group_names.index(None))
    return spc_group

def set_reactant_defaults(reactants):
    for reactant in reactants:
        reactant['sign'] = '-'
        if reactant['stoic'] is None:
            reactant['stoic'] = '1'
    return consolidate_dups(remove_nones(reactants))

def set_product_defaults(products):
    for product in products:
        if product['sign'] is None:
            product['sign'] = '+'
        if product['stoic'] is None:
            product['stoic'] = '1'
    product_names = [p['name'] for p in products]
    return consolidate_dups(remove_nones(products))

p = re.compile('(\d) ([A-Z])')

for ri,rxn in enumerate(reaction_list):
    reactants = p.sub(r'\1*\2',rxn["reactants"])

    products = p.sub(r'\1*\2',str(rxn["products"]).replace('None','').replace('\n',''))
    
    if products is None:
        products = ''
    
    if rxn['reaction_type'] == 'Photolysis':
        reaction_type = 'j'
    else:
        reaction_type = 'k'
    
    reaction_yaml+=(1*tab)+("RXN_%02d: %s =%s> %s" % (int(rxn["r_num"]),reactants,reaction_type,products)).replace('\n','').replace('  ',' ')+'\n'
        
    

        
    reactants = set_reactant_defaults( \
                  [i.groupdict() for i in \
                       species_pattern.finditer(reactants) \
                  ])

    products =  set_product_defaults( \
                  [i.groupdict() for i in \
                       species_pattern.finditer(products)
                  ])

    for role,spc_group in [("r", reactants), ("p", products)]: 
        for spc in spc_group: 
            spc["name"] = spc["name"].replace("\n","")
            if spc["name"] not in exclude_species:
                if spc["name"] not in species_list:
                    species_list.append(spc["name"])
                    species_yaml += 1 * tab + "- &%s '%s'\n" % (spc["name"],spc["name"])


initial_yaml = '---\n' + \
               species_yaml + '\n' +  \
               reaction_yaml + '\n' + \
               file('../../mechanisms/cb05_camx/cb05_camx_new_groups.yaml').read() + '\n...'

print >> file('../../mechanisms/cb05_camx/cb05_camx.yaml','wb'), initial_yaml
initial_mech = Mechanism(initial_yaml)
globals().update(initial_mech.species_dict)
net_reaction_rules=yaml.load(file('../../rules.yaml'))
net_reaction_yaml='\n'
net_reaction_yaml='net_reaction_list:\n'

def get_species(spc_str):
    return eval('['+spc_str+']',initial_mech.species_dict)
    
for net_key, net_rule in net_reaction_rules.iteritems():
    rxns = net_rule.get('rxns', [])
    if net_rule.has_key('r') or net_rule.has_key('p') or not net_rule.has_key('rxns'):
        reactants = get_species(net_rule.get('r',''))
        products = get_species(net_rule.get('p',''))
        logical_and = net_rule.get('logical_and',True)
        rxns += initial_mech.find_rxns(reactants,products,logical_and)
    
    
    net_reaction_yaml += '    %s: %s\n' % (net_key, ' + '.join(rxns))

full_yaml = '---\n' + \
               species_yaml + '\n' +  \
               reaction_yaml + '\n' + \
               net_reaction_yaml + '\n' + \
               file('../../mechanisms/cb05_camx/cb05_camx_new_groups.yaml').read() + '\n' + \
               file('camx_netphysical26.yaml').read() + '\n...'
               
# + '\n' + \
#               file('diagram.yaml').read()
print full_yaml