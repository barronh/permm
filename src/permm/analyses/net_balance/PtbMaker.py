def PtbTable(mech):
    o3 = mech('FCONC[O3].array()')
    peak_idx = (o3 == o3.max()).nonzero()
    template = '%-8s, %-15s\t%8.1f\t%8.1f\n'
    result = ''
    
    ## Nitrogen
    #Initial NOx
    spc = 'NOx'
    proc = 'Initial'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[0]
    peak = hourly[peak_idx]
    result += template % (spc, proc.lower(), daily, peak)
    
    #Avg NOx
    spc = 'NOx'
    proc = 'average'
    hourly = mech('(Initial[%s].array()+Final[%s].array())/2' % (spc,spc))
    daily = hourly.mean()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Avg NOy
    spc = 'NOy'
    proc = 'average'
    hourly = mech('Initial[%s].array()' % (spc,))
    daily = hourly.mean()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Emis NOx
    spc = 'NOx'
    proc = 'Emissions'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #H_Trans NOy
    spc = 'NOy'
    proc = 'H_Trans'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #H_Trans NOy
    spc = 'NOy'
    proc = 'V_Trans'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Entrain NOy
    spc = 'NOy'
    proc = 'Entrain'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Deposit NOy
    spc = 'NOy'
    proc = 'Deposit'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Final NOx
    spc = 'NOx'
    proc = 'Final'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[-1]
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Final NOy
    spc = 'NOy'
    proc = 'Final'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[-1]
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    ## Carbon
    #Initial VOCm
    spc = 'VOCm'
    proc = 'Initial'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[0]
    peak = hourly[peak_idx]
    result += template % (spc, proc.lower(), daily, peak)
    
    #Avg VOCm
    spc = 'VOCm'
    proc = 'average'
    hourly = mech('(Initial[%s].array()+Final[%s].array())/2' % (spc,spc))
    daily = hourly.mean()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Avg VOCC
    spc = 'VOCC'
    proc = 'average'
    hourly = mech('Initial[%s].array()' % (spc,))
    daily = hourly.mean()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Emis VOCm
    spc = 'VOCm'
    proc = 'Emissions'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #H_Trans VOCm
    spc = 'VOCm'
    proc = 'H_Trans'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #H_Trans VOCm
    spc = 'VOCm'
    proc = 'V_Trans'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Entrain VOCm
    spc = 'VOCm'
    proc = 'Entrain'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Deposit VOCm
    spc = 'VOCm'
    proc = 'Deposit'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly.sum()
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Final VOCm
    spc = 'VOCm'
    proc = 'Final'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[-1]
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    #Final VOCC
    spc = 'VOCC'
    proc = 'Final'
    hourly = mech('%s[%s].array()' % (proc, spc))
    daily = hourly[-1]
    peak = hourly[peak_idx]
    result += template % (spc, proc, daily, peak)
    
    ## VOC / NOx
    #Emissions VOCm/NOx
    numerator_spc = 'VOCm'
    denominator_spc = 'NOx'
    proc = 'Emissions'
    numerator_hourly = mech('%s[%s].array()' % (proc, numerator_spc))
    numerator_daily = numerator_hourly[-1]
    numerator_peak = numerator_hourly[peak_idx]
    denominator_hourly = mech('%s[%s].array()' % (proc, denominator_spc))
    denominator_daily = denominator_hourly[-1]
    denominator_peak = denominator_hourly[peak_idx]
    daily = numerator_daily.sum()/denominator_daily.sum()
    peak = numerator_peak/denominator_peak
    result += template % ('VOCm/NOx', proc, daily, peak)
    
    #Emissions VOCm/NOx
    numerator_spc = 'VOCm'
    denominator_spc = 'NOx'
    proc = 'average'
    numerator_hourly = mech('(Initial[%s].array()+Final[%s].array())/2' % (numerator_spc, numerator_spc))
    numerator_daily = numerator_hourly[-1]
    numerator_peak = numerator_hourly[peak_idx]
    denominator_hourly = mech('(Initial[%s].array()+Final[%s].array())/2' % (denominator_spc, denominator_spc))
    denominator_daily = denominator_hourly[-1]
    denominator_peak = denominator_hourly[peak_idx]
    daily = numerator_daily.sum()/denominator_daily.sum()
    peak = numerator_peak/denominator_peak
    result += template % ('VOCm/NOx', proc, daily, peak)
    
    return result
