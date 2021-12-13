from numpy import exp,sin,radians

def clrexc_multiplier(clrexc_dist,clrexc_b,clrexc_filt,clrexc_expcor):
    if (clrexc_filt=='J'):
        # standard, (B-V)_0 = (B-V) - E(B-V)
        return (1.0-clrexc_expcor*exp(-clrexc_dist*sin(abs(radians(clrexc_b)))/200.0))
    elif (clrexc_filt=='G'):
        #return (1.0-clrexc_expcor*exp(-clrexc_dist*sin(abs(radians(clrexc_b)))/200.0))*(1.072-0.634)/0.324
        return (1.0-clrexc_expcor*exp(-clrexc_dist*sin(abs(radians(clrexc_b)))/200.0))*1.60
    elif (clrexc_filt=='2'):
        return (1.0-clrexc_expcor*exp(-clrexc_dist*sin(abs(radians(clrexc_b)))/125.0))*(0.819-0.350)
    else:
        # used for debugging
        return (0.0)