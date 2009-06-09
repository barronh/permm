__all__ = ['SumTable']
from numpy import ndarray, where, nan

def sum_line(label, hourly, aggregate = None):
    if aggregate is None:
        aggregate = hourly.sum()
    nhrs = hourly.shape[0]
    hourly = hourly.tolist()
    
    template = '%-21s' + '%11.4f' * (nhrs + 1) + '\n'
    values = tuple([label] + hourly + [aggregate])
    
    return template % values

def time_line(hourly):
    nhrs = hourly.shape[0]
    hourly = hourly.tolist()
    
    template = '%-21s' % 'Time' + '% 11i' * nhrs + '%11s\n' % 'Sum'
    values = tuple(hourly)
    
    return template % values
    
def SumTable(mech):
    nrxns = mech.nreaction_dict
    globals().update(mech.species_dict)
    time = mech.mrg.variables['TFLAG'][:,0,1]
    def header(title):
        result = '\n%s\n' % (title,)
        result += time_line(time)
        return result
        
    result = ''
    inorg_pl_hv = nrxns['Inorganic New Radical']
    ald_pl_hv = nrxns['ORG+hv radical source']
    ox_pl_org = nrxns['Ox+org radical source']
    organic_new_rads = ald_pl_hv + ox_pl_org
    oh_pl_org = nrxns['OH+(organic+NO2)']
    total_new_rads = organic_new_rads + inorg_pl_hv
    total_prod_rads = oh_pl_org + total_new_rads
    pan_rads = nrxns['PAN Production']
    
    ### CXO3 New, Prod, and Loss
    cxo3_from_pan = pan_rads[CXO3]
    cxo3_non_pan = nrxns['CXO3 -> Radical'][CXO3]
    no2_from_cxo3 = nrxns['CXO3 -> Radical'][NO2]
    xo2_from_cxo3 = nrxns['CXO3 -> Radical'][XO2]
    ho2_from_cxo3 = nrxns['CXO3 -> Radical'][pHO2]
    no2_per_cxo3 = no2_from_cxo3 / total_prod_rads[CXO3]
    no2_per_cxo3_agg = no2_from_cxo3.sum() / total_prod_rads[CXO3].sum()
    xo2_per_cxo3 = xo2_from_cxo3 / total_prod_rads[CXO3]
    xo2_per_cxo3_agg = xo2_from_cxo3.sum() / total_prod_rads[CXO3].sum()
    ho2_per_cxo3 = ho2_from_cxo3 / total_prod_rads[CXO3]
    ho2_per_cxo3_agg = ho2_from_cxo3.sum() / total_prod_rads[CXO3].sum()
    prod_cxo3 = total_prod_rads[CXO3]
    
    result += header('CXO3 New, Prod, and Loss')
    # No photolysis production
    #result += sum_line('  Ald+hv', ald_pl_hv[CXO3])
    result += sum_line('  Ox+org (inc: NO3)', ox_pl_org[CXO3])
    result += sum_line('Total new CXO3', total_new_rads[CXO3])
    result += sum_line('OH+org', oh_pl_org[CXO3])
    result += sum_line('Total CXO3 Prod', prod_cxo3)
    result += sum_line('CXO3+NO2 PAN Prod', cxo3_from_pan)
    result += sum_line('CXO3+NO Loss', cxo3_non_pan)
    result += sum_line('  NO2', no2_from_cxo3)
    result += sum_line('  XO2', xo2_from_cxo3)
    result += sum_line('  HO2', ho2_from_cxo3)
    result += sum_line('Total CXO3 Loss', cxo3_from_pan + cxo3_non_pan)    
    result += sum_line('NO2/CXO3', no2_per_cxo3, no2_per_cxo3_agg)
    result += sum_line('XO2/CXO3', xo2_per_cxo3, xo2_per_cxo3_agg)
    result += sum_line('HO2/CXO3', ho2_per_cxo3, ho2_per_cxo3_agg)

    ### C2O3 New, Prod, and Loss
    c2o3_from_pan = pan_rads[C2O3]
    c2o3_non_pan = nrxns['C2O3 -> Radical'][C2O3]
    no2_from_c2o3 = nrxns['C2O3 -> Radical'][NO2]
    meo2_from_c2o3 = nrxns['C2O3 -> Radical'][MEO2]
    no2_per_c2o3 = no2_from_c2o3 / total_prod_rads[C2O3]
    no2_per_c2o3_agg = no2_from_c2o3.sum() / total_prod_rads[C2O3].sum()
    meo2_per_c2o3 = meo2_from_c2o3 / total_prod_rads[C2O3]
    meo2_per_c2o3_agg = meo2_from_c2o3.sum() / total_prod_rads[C2O3].sum()
    prod_c2o3 = total_prod_rads[C2O3]
    
    result += header('C2O3 New, Prod, and Loss')
    result += sum_line('  Ald+hv', ald_pl_hv[C2O3])
    result += sum_line('  Ox+org (inc: NO3)', ox_pl_org[C2O3])
    result += sum_line('Total new C2O3', total_new_rads[C2O3])
    result += sum_line('OH+org', oh_pl_org[C2O3])
    result += sum_line('Total C2O3 Prod', prod_c2o3)
    result += sum_line('C2O3+NO2 PAN Prod', c2o3_from_pan)
    result += sum_line('C2O3+NO Loss', c2o3_non_pan)
    result += sum_line('  NO2', no2_from_c2o3)
    result += sum_line('  MEO2', meo2_from_c2o3)
    result += sum_line('  -OOX', nrxns['C2O3 -> Radical'][RO2])
    result += sum_line('Total C2O3 Loss', c2o3_from_pan + c2o3_non_pan)
    result += sum_line('NO2/C2O3', no2_per_c2o3, no2_per_c2o3_agg)
    result += sum_line('MEO2/C2O3', meo2_per_c2o3, meo2_per_c2o3_agg)

    ### MEO2 New, Prod, and Loss
    prod_meo2 = meo2_from_c2o3+total_prod_rads[MEO2]
    ho2_from_meo2 = nrxns['MEO2 -> Radical'][pHO2]
    no2_from_meo2 = nrxns['MEO2 -> Radical'][NO2]
    no2_per_meo2 = no2_from_meo2 / prod_meo2
    no2_per_meo2_agg = no2_from_meo2.sum() / prod_meo2.sum()
    ho2_per_meo2 = ho2_from_meo2 / prod_meo2
    ho2_per_meo2_agg = ho2_from_meo2.sum() / prod_meo2.sum()

    result += header('MEO2 New, Prod, and Loss')
    result += sum_line('  Ald+hv', ald_pl_hv[MEO2])
    result += sum_line('Total new MEO2', total_new_rads[MEO2])
    result += sum_line('MEO2 from C2O3', meo2_from_c2o3)
    result += sum_line('OH+org', oh_pl_org[MEO2])
    result += sum_line('Total MEO2 Prod', prod_meo2)
    result += sum_line('MEO2+NO Loss', nrxns['MEO2 -> Radical'][MEO2])
    result += sum_line('  HO2', ho2_from_meo2)
    result += sum_line('  NO2', no2_from_meo2)
    result += sum_line('NO2/MEO2', no2_per_meo2, no2_per_meo2_agg)
    result += sum_line('HO2/MEO2', ho2_per_meo2, ho2_per_meo2_agg)
    
    ### XO2 New, Prod, and Loss
    prod_xo2 = total_prod_rads[XO2_TO2] + xo2_from_cxo3
    prod_xo2n = total_prod_rads[XO2N]
    prod_xo2_n = prod_xo2 + prod_xo2n
    ho2_from_xo2 = nrxns['Whole_Mechanism'][xHO2]
    no2_per_xo2_n = hourly = nrxns['XO2_N_TO2 -> Radical'][NO2] / prod_xo2_n
    no2_per_xo2_n_agg = agg    = nrxns['XO2_N_TO2 -> Radical'][NO2].sum() / prod_xo2_n.sum()
    ho2_per_xo2_n = hourly = ho2_from_xo2 / prod_xo2_n
    ho2_per_xo2_n_agg =  ho2_from_xo2.sum() / prod_xo2_n.sum()

    result += header('XO2 New, Prod, and Loss')
    result += sum_line('  Ald+hv', ald_pl_hv[XO2_TO2])
    result += sum_line('  Ox+org (inc: NO3)', ox_pl_org[XO2_TO2])
    result += sum_line('Total new XO2', total_new_rads[XO2_TO2])
    result += sum_line('XO2 from CXO3', xo2_from_cxo3)
    result += sum_line('OH+org', oh_pl_org[XO2_TO2])
    result += sum_line('Total XO2 Prod', prod_xo2)
    result += sum_line('  Ox+org (inc: NO3)', ox_pl_org[XO2N])
    result += sum_line('Total new XO2N', organic_new_rads[XO2N])
    result += sum_line('OH+org', oh_pl_org[XO2N])
    result += sum_line('Total XO2N Prod', prod_xo2n)
    result += sum_line('Total XO2/XO2N Prod', prod_xo2_n)
    result += sum_line('  XO2+NO Loss', nrxns['XO2_TO2 -> Radical'][XO2_TO2])
    result += sum_line('  XO2N+NO Loss', nrxns['XO2N -> Radical'][XO2N])
    result += sum_line('Total XO2/XO2N Loss', nrxns['XO2_N_TO2 -> Radical'][XO2_N_TO2])
    result += sum_line('  NO', nrxns['XO2_N_TO2 -> Radical'][NO])
    result += sum_line('  NO2', nrxns['XO2_N_TO2 -> Radical'][NO2])
    result += sum_line('  HO2', ho2_from_xo2)
    result += sum_line('  NTR', nrxns['XO2_N_TO2 -> Radical'][NTR])
    result += sum_line('  -OOX', nrxns['XO2_N_TO2 -> Radical'][RO2])
    result += sum_line('NO2/XO2', no2_per_xo2_n, no2_per_xo2_n_agg)
    result += sum_line('HO2/XO2', ho2_per_xo2_n, ho2_per_xo2_n_agg)

    ### HO2 New, Prod, and Loss
    ho2_from_hco3 = nrxns['HCO3 -> Radical'][HO2]
    prod_ho2 = total_prod_rads[pHO2] + ho2_from_meo2 + ho2_from_xo2 + ho2_from_hco3
    no2_from_ho2 = nrxns['HO2 -> Radical'][NO2]
    oh_from_ho2 = nrxns['HO2 -> Radical'][OH]
    no2_per_ho2 = no2_from_ho2 / prod_ho2
    no2_per_ho2_agg = no2_from_ho2.sum() / prod_ho2.sum()
    oh_per_ho2 = oh_from_ho2 / prod_ho2
    oh_per_ho2_agg = oh_from_ho2.sum() / prod_ho2.sum()

    result += header('HO2 New, Prod, and Loss')
    result += sum_line('Inorg+hv', inorg_pl_hv[pHO2])
    result += sum_line('  Ald+hv', ald_pl_hv[pHO2])
    result += sum_line('  Ox+org (inc: NO3)', ox_pl_org[pHO2])
    result += sum_line('Total new HO2', total_new_rads[pHO2])
    result += sum_line('HO2 from MEO2', ho2_from_meo2)
    result += sum_line('HO2 from XO2', ho2_from_xo2)
    result += sum_line('HO2 from HCO3', ho2_from_hco3)
    result += sum_line('OH+org', oh_pl_org[pHO2])
    result += sum_line('Total HO2 Prod', prod_ho2)
    result += sum_line('  CXO3', nrxns['RXN_108'][lHO2])
    result += sum_line('  O3', nrxns['RXN_13'][lHO2])
    result += sum_line('  CRO', nrxns['RXN_137'][lHO2])
    result += sum_line('  HO2', nrxns['RXN_34'][lHO2])
    result += sum_line('  H2O', nrxns['RXN_35'][lHO2])
    result += sum_line('  OH', nrxns['RXN_43'][lHO2])
    result += sum_line('  O', nrxns['RXN_44'][lHO2])
    result += sum_line('  NO3', nrxns['RXN_48'][lHO2])
    result += sum_line('  XO2/XO2N', (nrxns['RXN_56']+nrxns['RXN_57'])[lHO2])
    result += sum_line('  MEO2', nrxns['RXN_69'][lHO2])
#    result += sum_line('  FORM', nrxns['RXN_79'][lHO2])
#    result += sum_line('  HCO3', nrxns['RXN_82'][lHO2])
    result += sum_line('  C2O3', nrxns['RXN_92'][lHO2])
    result += sum_line('  NO', nrxns['RXN_30'][lHO2])
    result += sum_line('    NO2', no2_from_ho2)
    result += sum_line('    OH', oh_from_ho2)
    result += sum_line('Total HO2 Loss', nrxns['Whole_Mechanism'][lHO2])
    result += sum_line('NO2/HO2',no2_per_ho2,no2_per_ho2_agg)
    result += sum_line('OH/HO2',oh_per_ho2,oh_per_ho2_agg)

    ### OH New, Prod, and Loss
    oh_per_cxo3 = ( \
                    xo2_per_cxo3 * ho2_per_xo2_n + \
                    ho2_per_cxo3 \
                   ) * \
                   oh_per_ho2
    oh_from_new_cxo3 = total_new_rads[CXO3] * oh_per_cxo3
    oh_from_reused_cxo3 = (prod_cxo3 - total_new_rads[CXO3]) * oh_per_cxo3
                   
    oh_per_c2o3 = meo2_per_c2o3 * ho2_per_meo2 * oh_per_ho2
    oh_from_new_c2o3 = total_new_rads[C2O3] * oh_per_c2o3
    oh_from_reused_c2o3 = (prod_c2o3 - total_new_rads[C2O3]) * oh_per_c2o3

    oh_per_meo2 = ho2_per_meo2 * oh_per_ho2
    oh_from_new_meo2 = total_new_rads[MEO2] * oh_per_meo2
    oh_from_reused_meo2 = (prod_meo2 - total_new_rads[MEO2]) * oh_per_meo2
    
    oh_per_xo2 = ho2_per_xo2_n * oh_per_ho2
    oh_from_new_xo2_n = total_new_rads[XO2_N_TO2] * oh_per_xo2
    oh_from_reused_xo2_n = (prod_xo2_n - total_new_rads[XO2_N_TO2]) * oh_per_xo2

    oh_from_new_ho2 = total_new_rads[pHO2] * oh_per_ho2
    oh_from_reused_ho2 = (prod_ho2 - total_new_rads[pHO2]) * oh_per_ho2

    oh_from_new_rad = oh_from_new_cxo3 + oh_from_new_c2o3 + oh_from_new_xo2_n + \
                      oh_from_new_meo2 + oh_from_new_ho2
    oh_new_inorg = inorg_pl_hv[OH]
    
    oh_total_new = oh_from_new_rad + oh_new_inorg
    
    oh_total = nrxns['Total OH Reacted'][OH]
    oh_from_reused_rad = -oh_total - oh_total_new
    
    oh_chain_length = -oh_total / oh_total_new
    oh_chain_length_agg = -oh_total.sum() / oh_total_new.sum()
    
    oh_propef = 1 - 1/oh_chain_length
    oh_propef_agg = 1 - 1/oh_chain_length_agg

    oh_reacted_with_co = -mech.make_net_rxn([OH,CO],[])[OH]

    oh_reacted_with_voc = -mech.make_net_rxn([OH,VOC],[])[OH]
    
    oh_reacted_with_hc = -mech.make_net_rxn([OH, HC],[])[OH]

    oh_reacted_with_no2 = -mech.make_net_rxn([OH,NO2],[])[OH]

    result += header('OH New, Prod, and Loss')
    result += sum_line('  nOH from CXO3', oh_from_new_cxo3)
    result += sum_line('  nOH from C2O3', oh_from_new_c2o3)
    result += sum_line('  nOH from XO2', oh_from_new_xo2_n)
    result += sum_line('  nOH from MEO2', oh_from_new_meo2)
    result += sum_line('  nOH from HO2', oh_from_new_ho2)
    result += sum_line('nOH from new Rad', oh_from_new_rad)
    result += sum_line('nOH from Inorg+hv', oh_new_inorg)
    result += sum_line('nOH Total', oh_total_new)
    result += sum_line('OH reused Rad (diff)', oh_from_reused_rad)
    result += sum_line('OH Total Reacted', oh_total)
    result += sum_line('OH Chain length', oh_chain_length, oh_chain_length_agg)
    result += sum_line('OH Prop Eff', oh_propef, oh_propef_agg)
    result += sum_line('OH + CO', oh_reacted_with_co)
    result += sum_line('OH + VOC', oh_reacted_with_voc)
    result += sum_line('OH + HC', oh_reacted_with_hc)
    result += sum_line('OH + NO2', oh_reacted_with_no2)

    ### NO2 New, Prod, and Loss
    no2_emissions = mech('Emissions[NO2].array()')
    no2_htrans = mech('H_Trans[NO2].array()')
    no2_htrans_gain = where(no2_htrans > 0, no2_htrans, 0)
    no2_htrans_loss = where(no2_htrans < 0, no2_htrans, 0)
    no2_vtrans = mech('V_Trans[NO2].array()')
    no2_vtrans_gain = where(no2_vtrans > 0, no2_vtrans, 0)
    no2_vtrans_loss = where(no2_vtrans < 0, no2_vtrans, 0)
    no2_motion = mech('Motion[NO2].array()')
    no2_motion_gain = where(no2_motion > 0, no2_motion, 0)
    no2_motion_loss = where(no2_motion < 0, no2_motion, 0)
  
    no2_initial = mech('Initial[NO2]').array()
    no2_carry  = mech('Initial[NO2]').array()
    no2_initial[1:] = 0
    no2_carry[:1] = 0

    no2_new = no2_initial + no2_emissions + no2_htrans_gain + \
              no2_vtrans_gain + no2_motion_gain

    result += header('New NO2, Srcs Outside Box')
    result += sum_line('  NO2 initial', no2_initial)
    result += sum_line('  NO2 emissions', no2_emissions)
    result += sum_line('  NO2 horiz trans', no2_htrans_gain)
    result += sum_line('  NO2 vert trans', no2_vtrans_gain)
    result += sum_line('  NO2 entrain/dil', no2_motion_gain)
    result += sum_line('New NO2', no2_new)

    no2_phot = nrxns['NO2 Photolysis1']
    no2_from_o3_oxid = where(no2_phot[O3]<0, no2_phot[NO2], 0)
    no2_from_ro2_oxid = nrxns['NO+RO2 Oxidation'][NO2]
    no2_from_no3_pl_org = nrxns['NO3+org radical source'][NO2]
    
    no2_available = no2_initial + no2_carry + no2_from_o3_oxid + no2_new + no2_from_ro2_oxid


    # Available[NO2] - CHEM[NO2] - PHY[NO2] = 0
    # only if no2_available includes no2_initial
    no2_available_agg = (no2_initial + no2_from_o3_oxid + no2_new + no2_from_ro2_oxid).sum()
    
    no2_from_chem = no2_from_ro2_oxid + no2_from_no3_pl_org
    no2_frac_from_chem = no2_from_chem / no2_available
    
    no2_frac_from_chem[0] = no2_from_chem[0] / (no2_available[0] - no2_initial[0])

    # Skipping assumptions
    # Assuming initial NO2 chemicaly produced fraction is 
    # proportional to first hour process fractions
    #no2_from_chem[0] += no2_frac_from_chem[0] * no2_initial[0]

    no2_frac_from_chem_agg = no2_from_chem.sum() / no2_available_agg

    for hr in range(no2_from_chem.shape[0])[1:]:
        no2_from_chem[hr] += no2_frac_from_chem[hr-1] * no2_carry[hr]
        no2_frac_from_chem = no2_from_chem / no2_available
    result += header('Total Available NO2 Inside Box')
    result += sum_line('Initial', no2_initial)
    result += sum_line('Carry over', no2_carry, 0)
    result += sum_line('O3+NO Oxid', no2_from_o3_oxid)
    result += sum_line('New NO2', no2_new)
    result += sum_line('(RO2+Ox)+NO Oxid',no2_from_ro2_oxid + no2_from_no3_pl_org)
    result += sum_line('Total Available NO2', no2_available, no2_available_agg)
    result += sum_line('Frac NO2 from chem', no2_frac_from_chem, no2_frac_from_chem_agg)
    
    no2_final = -mech('Final[NO2]').array()
    no2_dep = mech('Deposit[NO2]').array()
    
    no2_phy_loss = no2_final + no2_dep + no2_htrans_loss + no2_vtrans_loss + no2_motion_loss
    no2_phy_loss_agg = no2_final[-1] + (no2_dep + no2_htrans_loss + no2_vtrans_loss + no2_motion_loss).sum()
    result += header('Physical NO2 Losses')
    result += sum_line('  NO2 Final', no2_final, no2_final[-1])
    result += sum_line('  NO2 deposition', no2_dep)
    result += sum_line('  NO2 horiz trans', no2_htrans_loss)
    result += sum_line('  NO2 vert  trans', no2_vtrans_loss)
    result += sum_line('  NO2 entrain/dil', no2_motion_loss)
    result += sum_line('Total NO2 phy losses', no2_phy_loss, no2_phy_loss_agg)
    result += sum_line('Total NO2 for Chem', no2_phy_loss + no2_available, no2_phy_loss_agg + no2_available_agg)
    
    no2_phot_loss = where(no2_phot[NO2]<0, no2_phot[NO2], 0)
    no_prod = where(no2_phot[NO2]<0, no2_phot[NO], 0)
    o3_prod = where(no2_phot[NO2]<0, no2_phot[O3], 0)
    no3_prod = where(no2_phot[NO2]<0, no2_phot[NO3], 0)
    
    o3_per_no2 = o3_prod / no2_available
    o3_per_no2_agg = o3_prod.sum() / no2_available_agg

    no_per_no2 = no_prod / no2_available
    no_per_no2_agg = no_prod.sum() / no2_available_agg
    
    no2_to_pans = nrxns['PAN Production'][NO2]
    no2_to_ntr = -nrxns['NO2 Termination'][NTR]
    no2_to_hno3 = -nrxns['NO2 Termination'][HNO3]
    no2_term = no2_to_pans + no2_to_ntr + no2_to_hno3
    no_new_sec = no2_new * no_per_no2

    no_initial = mech('Initial[NO]').array()
    no_carry  = mech('Initial[NO]').array()
    no_initial[1:] = 0
    no_carry[:1] = 0
    no_emissions = mech('Emissions[NO].array()')
    no_htrans = mech('H_Trans[NO].array()')
    no_htrans_gain = where(no_htrans > 0, no_htrans, 0)
    no_htrans_loss = where(no_htrans < 0, no_htrans, 0)
    no_vtrans = mech('V_Trans[NO].array()')
    no_vtrans_gain = where(no_vtrans > 0, no_vtrans, 0)
    no_vtrans_loss = where(no_vtrans < 0, no_vtrans, 0)
    no_motion = mech('Motion[NO].array()')
    no_motion_gain = where(no_motion > 0, no_motion, 0)
    no_motion_loss = where(no_motion < 0, no_motion, 0)
    no_from_o3_oxid = where(no2_phot[O3]<0, no2_phot[NO], 0)
    
    no_total_new = no_initial + no_carry + no_from_o3_oxid + \
                   no_emissions + no_htrans_gain + \
                   no_vtrans_gain + no_motion_gain + \
                   no_new_sec

    no_total_new_agg = (no_initial +  no_from_o3_oxid + \
                   no_emissions + no_htrans_gain + \
                   no_vtrans_gain + no_motion_gain + \
                   no_new_sec).sum()

    no_chain_length = no_prod / no_total_new
    no_chain_length_agg = no_prod.sum() / no_total_new_agg

    no_recreated = no_prod - no_new_sec
    no_recreated_agg = no_prod.sum() - no_new_sec.sum()
    no_propef = 1 - 1/where(no_chain_length == 0, nan, no_chain_length)
    no_propef_agg = 1 - 1/no_chain_length_agg


    result += header('Chemical NO2 Losses')
    result += sum_line('   HNO3 formation', no2_to_hno3)
    result += sum_line('   PAN  formation', no2_to_pans)
    result += sum_line('   NTR  formation', no2_to_ntr)
    result += sum_line('Total NO2 Term Loss', no2_term)
    result += sum_line('   NO2 Photolyzed', no2_phot_loss)
    result += sum_line('   NO  Produced', no_prod)
    result += sum_line('   O3  Produced', o3_prod)
    result += sum_line('   NO3 Produced', no3_prod)
    result += sum_line('Total NO2 Chem Loss', no2_term + no2_phot_loss)
    result += sum_line('O3 yld per NO2 phot', o3_per_no2, o3_per_no2_agg)
    result += sum_line('NO yld per NO2 avail', no_per_no2, no_per_no2_agg)
    result += sum_line('new NO frm new NO2', no_new_sec)
    result += sum_line('recreated NO', no_recreated, no_recreated_agg)

    ### New NO, Physical NO sources

    result += header('New NO, Physical NO sources')  
    result += sum_line('  Initial', no_initial)
    result += sum_line('  Carry over', no_carry, 0)
    result += sum_line('  O3+NO Oxid', no_from_o3_oxid)
    result += sum_line('  NO emissions', no_emissions)
    result += sum_line('  NO horiz trans', no_htrans_gain)
    result += sum_line('  NO vert trans', no_vtrans_gain)
    result += sum_line('  NO entrain/dil', no_motion_gain)
    result += sum_line('  sec new NO', no_new_sec)
    result += sum_line('Total new NO', no_total_new, no_total_new_agg)
    result += sum_line('NO P_n', no_propef, no_propef_agg)
    result += sum_line('NO chain length', no_chain_length, no_chain_length_agg)
    
    return result
