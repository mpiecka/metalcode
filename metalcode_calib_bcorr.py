from metalcode_v1_2_0_calib_tempe import Teff

def BolCorBV(BC_BV,BC_FeH,BC_label):
    BolCor=-999.9
    if (BC_label=='J'):
        #if ((BC_BV<0.1) and (BC_BV>-0.2)):
        #    BolCor=-0.155+2.84*BC_BV-15.13*BC_BV**2+25.57*BC_BV**3
        #elif ((BC_BV<0.8) and (BC_BV>0.1)):
        #    BolCor=-0.030+0.491*BC_BV-0.937*BC_BV**2+0.055*BC_BV**3
        if ((BC_BV<0.8) and (BC_BV>0.2)):
            BolCor=-0.030+0.491*BC_BV-0.937*BC_BV**2+0.055*BC_BV**3
    elif (BC_label=='G'):
        BC_Teff=Teff(BC_BV,BC_BV,BC_FeH,BC_label)
        BC_coeff1=[6.000e-2,6.731e-5,-6.647e-8,2.859e-11,-7.197e-15]
        BC_coeff2=[1.749e0,1.977e-3,3.737e-7,-8.966e-11,-4.183e-14]
        if ((BC_Teff<8000.0) and (BC_Teff>4000.0)):
            BolCor=0.0
            for BC_i in range(len(BC_coeff1)):
                BolCor+=BC_coeff1[BC_i]*(BC_Teff - 5772.0)**(BC_i)
        elif ((BC_Teff<4000.0) and (BC_Teff>3300.0)):
            BolCor=0.0
            for BC_i in range(len(BC_coeff2)):
                BolCor+=BC_coeff2[BC_i]*(BC_Teff - 5772.0)**(BC_i)
    elif (BC_label=='2'):
        if ((BC_BV<0.80) and (BC_BV>0.10)):
            BolCor=0.366+2.485*BC_BV-1.232*BC_BV**2
    return BolCor

"""
def BolCorTeff(BC_Teff):
    BolCor=-999.9
    if ((log10(BC_Teff)<4.25) and (log10(BC_Teff)>4.00)):
        BolCor=-0.238-5.561*(log10(BC_Teff)-4.0)
    elif ((log10(BC_Teff)<4.00) and (log10(BC_Teff)>3.76)):
        BolCor=-0.243-5.0*(log10(BC_Teff)-4.0)-26.0*(log10(BC_Teff)-4.0)**2-33.33*(log10(BC_Teff)-4.0)**3 
    return BolCor
"""