from numpy import log10,exp,sin,radians

def absmag(absmag_E,absmag_m,absmag_d,absmag_label,absmag_b,absmag_expcor):
    #AxAy=1.324
    if (absmag_label=='J'):
        return (absmag_m-5.0*(log10(absmag_d)-1.0)-(1.0-absmag_expcor*exp(-absmag_d*sin(abs(radians(absmag_b)))/200.0))*absmag_E/0.324)
    elif (absmag_label=='G'):
        return (absmag_m-5.0*(log10(absmag_d)-1.0)-(1.0-absmag_expcor*exp(-absmag_d*sin(abs(radians(absmag_b)))/200.0))*0.832*absmag_E/0.324)
    elif (absmag_label=='2'):
        return (absmag_m-5.0*(log10(absmag_d)-1.0)-(1.0-absmag_expcor*exp(-absmag_d*sin(abs(radians(absmag_b)))/125.0))*0.819*absmag_E)