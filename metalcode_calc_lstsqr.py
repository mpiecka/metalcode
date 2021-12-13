from numpy import sqrt
from numpy import polyfit
from numpy import median
from numpy import std

def LstSqr(Lst_y,Lst_x,Lst_w,Lst_L,Lst_T): 
    # used for artificial weights, based on the dispersion in the LTN diagram  
    Lst_p=polyfit(Lst_y,Lst_x,4)
    Lst_z=[]
    for Lst_i in range(len(Lst_x)):
        Lst_a=0.0
        for Lst_j in range(len(Lst_p)):
            Lst_a+=Lst_y[Lst_i]**Lst_j*Lst_p[len(Lst_p)-1-Lst_j]
        Lst_z.append(abs(Lst_a-Lst_x[Lst_i]))

    # artificial weights to different parts of the LTN diagram
    Lst_q=[]
    for Lst_i in range(len(Lst_x)):
        #"""
        if (Lst_y[Lst_i]>3.00):
            Lst_q.append(0.10)
        elif (Lst_y[Lst_i]>1.00):
            Lst_q.append((0.012/median(Lst_z))**2.0)
        elif (Lst_y[Lst_i]>0.50):
            Lst_q.append((0.012/median(Lst_z))**1.4)
        elif (Lst_y[Lst_i]>-1.00):
            Lst_q.append(2.00)
        else:
            #Lst_q.append(0.0001)
            Lst_q.append(0.01)
        #"""
        #Lst_q.append(1.0)

    Lst_Yscale=0.1
    Lst_sum=[]
    for Lst_i in range(len(Lst_x)):
        Lst_rad=[]
        for Lst_j in range (len(Lst_T)):
            Lst_rad.append((Lst_x[Lst_i]-Lst_T[Lst_j])**2+(Lst_y[Lst_i]-Lst_L[Lst_j])**2)
        Lst_J=Lst_rad.index(min(Lst_rad))
        
        Lst_r0x=Lst_x[Lst_i]-Lst_T[Lst_J]
        Lst_r0y=Lst_y[Lst_i]-Lst_L[Lst_J]
        Lst_D0D=sqrt(Lst_r0x**2+Lst_Yscale*Lst_r0y**2)
        
        Lst_D1D=1.0e20
        if (Lst_J>0):
            Lst_r1x=Lst_T[Lst_J-1]-Lst_T[Lst_J]
            Lst_r1y=Lst_L[Lst_J-1]-Lst_L[Lst_J]
            Lst_d1x=Lst_T[Lst_J] + (Lst_r0x*Lst_r1x+Lst_r0y*Lst_r1y)/(1.0e-15+Lst_r1x**2+Lst_r1y**2)*Lst_r1x
            Lst_d1y=Lst_L[Lst_J] + (Lst_r0x*Lst_r1x+Lst_r0y*Lst_r1y)/(1.0e-15+Lst_r1x**2+Lst_r1y**2)*Lst_r1y
            if (sqrt(Lst_r1x**2+Lst_Yscale*Lst_r1y**2) >= sqrt((Lst_T[Lst_J]-Lst_d1x)**2+Lst_Yscale*(Lst_L[Lst_J]-Lst_d1y)**2)+sqrt((Lst_T[Lst_J-1]-Lst_d1x)**2+Lst_Yscale*(Lst_L[Lst_J-1]-Lst_d1y)**2)):
                Lst_D1D=sqrt((Lst_x[Lst_i]-Lst_d1x)**2+Lst_Yscale*(Lst_y[Lst_i]-Lst_d1y)**2)
        
        Lst_D2D=1.0e20  
        if (Lst_J<len(Lst_T)-1):  
            Lst_r2x=Lst_T[Lst_J+1]-Lst_T[Lst_J]
            Lst_r2y=Lst_L[Lst_J+1]-Lst_L[Lst_J]       
            Lst_d2x=Lst_T[Lst_J] + (Lst_r0x*Lst_r2x+Lst_r0y*Lst_r2y)/(1.0e-15+Lst_r2x**2+Lst_r2y**2)*Lst_r2x
            Lst_d2y=Lst_L[Lst_J] + (Lst_r0x*Lst_r2x+Lst_r0y*Lst_r2y)/(1.0e-15+Lst_r2x**2+Lst_r2y**2)*Lst_r2y
            if (sqrt(Lst_r2x**2+Lst_Yscale*Lst_r2y**2) >= sqrt((Lst_T[Lst_J]-Lst_d2x)**2+Lst_Yscale*(Lst_L[Lst_J]-Lst_d2y)**2)+sqrt((Lst_T[Lst_J+1]-Lst_d2x)**2+Lst_Yscale*(Lst_L[Lst_J+1]-Lst_d2y)**2)):
                Lst_D2D=sqrt((Lst_x[Lst_i]-Lst_d2x)**2+Lst_Yscale*(Lst_y[Lst_i]-Lst_d2y)**2)

        Lst_sum.append(min([Lst_D0D,Lst_D1D,Lst_D2D]))
        
    if (len(Lst_sum)>5):
        weightSum=0.0
        weightVar=0.0
        weightDen=0.0
        for Lst_i in range(len(Lst_sum)):
            weightSum+=Lst_w[Lst_i]*Lst_q[Lst_i]*Lst_sum[Lst_i]
            weightDen+=Lst_w[Lst_i]*Lst_q[Lst_i]
        weightSum=weightSum/float(weightDen)
        for Lst_i in range(len(Lst_sum)):
            weightVar+=Lst_w[Lst_i]*Lst_q[Lst_i] * (Lst_sum[Lst_i]-weightSum)**2
        weightVar=sqrt(weightVar/float(weightDen-1))
        return sqrt(weightSum*weightVar)
    else:
        return 1.0e20

def LstSca(Lst_y,Lst_x,Lst_L,Lst_T):
    Lst_Yscale=1.0
    Lst_sca=[]
    for Lst_i in range(len(Lst_x)):
        Lst_rad=[]
        for Lst_j in range (len(Lst_T)):
            Lst_rad.append((Lst_x[Lst_i]-Lst_T[Lst_j])**2+(Lst_y[Lst_i]-Lst_L[Lst_j])**2)
        Lst_J=Lst_rad.index(min(Lst_rad))
        
        Lst_r0x=Lst_x[Lst_i]-Lst_T[Lst_J]
        Lst_r0y=Lst_y[Lst_i]-Lst_L[Lst_J]
        Lst_D0D=sqrt(Lst_r0x**2+Lst_Yscale*Lst_r0y**2)
        Lst_D0S=Lst_r0x
        
        Lst_D1D=1.0e20
        Lst_D1S=1.0e20
        if (Lst_J>0):
            Lst_r1x=Lst_T[Lst_J-1]-Lst_T[Lst_J]
            Lst_r1y=Lst_L[Lst_J-1]-Lst_L[Lst_J]
            Lst_d1x=Lst_T[Lst_J] + (Lst_r0x*Lst_r1x+Lst_r0y*Lst_r1y)/(1.0e-15+Lst_r1x**2+Lst_r1y**2)*Lst_r1x
            Lst_d1y=Lst_L[Lst_J] + (Lst_r0x*Lst_r1x+Lst_r0y*Lst_r1y)/(1.0e-15+Lst_r1x**2+Lst_r1y**2)*Lst_r1y
            if (sqrt(Lst_r1x**2+Lst_Yscale*Lst_r1y**2) >= sqrt((Lst_T[Lst_J]-Lst_d1x)**2+Lst_Yscale*(Lst_L[Lst_J]-Lst_d1y)**2)+sqrt((Lst_T[Lst_J-1]-Lst_d1x)**2+Lst_Yscale*(Lst_L[Lst_J-1]-Lst_d1y)**2)):
                Lst_D1D=sqrt((Lst_x[Lst_i]-Lst_d1x)**2+Lst_Yscale*(Lst_y[Lst_i]-Lst_d1y)**2)
                Lst_D1S=(Lst_x[Lst_i]-Lst_d1x)
        
        Lst_D2D=1.0e20  
        Lst_D2S=1.0e20
        if (Lst_J<len(Lst_T)-1):  
            Lst_r2x=Lst_T[Lst_J+1]-Lst_T[Lst_J]
            Lst_r2y=Lst_L[Lst_J+1]-Lst_L[Lst_J]       
            Lst_d2x=Lst_T[Lst_J] + (Lst_r0x*Lst_r2x+Lst_r0y*Lst_r2y)/(1.0e-15+Lst_r2x**2+Lst_r2y**2)*Lst_r2x
            Lst_d2y=Lst_L[Lst_J] + (Lst_r0x*Lst_r2x+Lst_r0y*Lst_r2y)/(1.0e-15+Lst_r2x**2+Lst_r2y**2)*Lst_r2y
            if (sqrt(Lst_r2x**2+Lst_Yscale*Lst_r2y**2) >= sqrt((Lst_T[Lst_J]-Lst_d2x)**2+Lst_Yscale*(Lst_L[Lst_J]-Lst_d2y)**2)+sqrt((Lst_T[Lst_J+1]-Lst_d2x)**2+Lst_Yscale*(Lst_L[Lst_J+1]-Lst_d2y)**2)):
                Lst_D2D=sqrt((Lst_x[Lst_i]-Lst_d2x)**2+Lst_Yscale*(Lst_y[Lst_i]-Lst_d2y)**2)
                Lst_D2S=(Lst_x[Lst_i]-Lst_d2x)

        Lst_dss=[Lst_D0S,Lst_D1S,Lst_D2S]
        Lst_sum=[Lst_D0D,Lst_D1D,Lst_D2D]
        Lst_sca.append(Lst_dss[Lst_sum.index(min(Lst_sum))])

    return std(Lst_sca)