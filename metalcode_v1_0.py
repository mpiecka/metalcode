from numpy import *
import matplotlib.pyplot as plt
import time
import os
import warnings
warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)

# from metalcode_calib_metal import metal_transf
from metalcode_calib_tempe import Teff
from metalcode_calib_tempe import BolCorBV
from metalcode_calib_absmg import absmag
from metalcode_calib_clrex import clrexc_multiplier
from metalcode_calc_lstsqr import LstSqr




#-------------------------------------------------------
#---------------------INITIALISATION--------------------
#-------------------------------------------------------
def initialise():
    # input values - photometric system
    print('Initialisation ...')
    inputing=True
    list_vals=['G','J','2']
    while inputing:
        photosystem=input(' -- Pick photometric system (G,2,J): ')
        if (photosystem in list_vals):
            inputing=False
    # input values - grid spacing (age)
    inputing=True
    list_vals=[0.1,0.2]
    while inputing:
        try:
            age_step=float(input(' -- Pick grid spacing, age (0.1,0.2): '))
            if (age_step in list_vals):
                inputing=False
        except ValueError:
            pass
    # input values - grid spacing (Z)
    inputing=True
    list_vals=[0.005]
    while inputing:
        try:
            z_step=float(input(' -- Pick grid spacing, Z (0.005): '))
            if (z_step in list_vals):
                inputing=False
        except ValueError:
            pass
    # input values - Nredd
    inputing=True
    while inputing:
        try:
            Nredd=int(input(' -- Nredd (3): '))
            if ((Nredd == 0) or (Nredd % 2 == 1)):
                inputing=False
        except ValueError:
            pass
    # input values - reddening range
    inputing=True
    list_vals=[0.005]
    while inputing:
        try:
            redAdj=float(input(' -- Reddening range (0.0 .. 1.0): '))
            if (redAdj>0.0 and redAdj<1.0):
                inputing=False
        except ValueError:
            pass
    # input values - Niter
    inputing=True
    while inputing:
        try:
            Niter=int(input(' -- Niter (6): '))
            if (Niter >= 3):
                inputing=False
        except ValueError:
            pass
    return (photosystem,age_step,z_step,Nredd,Niter,redAdj)
#-------------------------------------------------------
#-------------------END INITIALISATION------------------
#-------------------------------------------------------





#-------------------------------------------------------
#-----------------------CALCULATIONS--------------------
#-------------------------------------------------------
# systematic corrections to the temperature calibration of observed data
#
# what you need is only q3 and q4 collections for all [Age,Z] in the grid,
# then you can correct for the systematics, works very well for less bright stars,
# but also slightly improves the situation for giants (especially for higher ages)
def sys_temp(sys_age,sys_z,sys_photosystem):
    global isochronesLTN
    global isochronesCMD
    global age_values
    global z_values

    sys_q1=[]
    sys_q2=[]
    sys_i=age_values.index(round(sys_age,1))
    sys_j=z_values.index(round(sys_z,3))
    sys_b,sys_a=isochronesLTN[sys_i][sys_j]
    sys_y,sys_x=isochronesCMD[sys_i][sys_j]
    sys_xx=[]
    sys_yy=[]
    sys_aa=[]
    sys_bb=[]
    for sys_k in range(len(sys_x)):
        sys_y0=logL( sys_y[sys_k] , BolCorBV(sys_x[sys_k],sys_y[sys_k],z_values[sys_j],sys_photosystem) )
        sys_x0=logTN( Teff(sys_x[sys_k],sys_y[sys_k],z_values[sys_j],sys_photosystem) , sys_y0 )
        if (sys_x0>-3 and sys_x0<3 and sys_y0>-5 and sys_y0<5):
            sys_yy.append(sys_y0)
            sys_xx.append(sys_x0)
            sys_aa.append(sys_a[sys_k])
            sys_bb.append(sys_b[sys_k])
            sys_q2.append(sys_a[sys_k]-sys_x0)
            sys_q1.append(sys_y0)

    sys_qx=sorted(sys_q1)
    sys_qy=[]
    for sys_i in range(len(sys_qx)):
        sys_qy.append(sys_q2[sys_q1.index(sys_qx[sys_i])])
    sys_q3=[]
    sys_q4=[]
    for sys_j in range(35):
        sys_qq=[]
        for sys_i in range(len(sys_qx)):
            if (sys_qx[sys_i]>(-1.0+sys_j*0.2) and sys_qx[sys_i]<=(-1.0+(sys_j+1)*0.2)):
                sys_qq.append(sys_qy[sys_i])
        if (len(sys_qq)>0):
            sys_q3.append((-1.0+(sys_j+0.5)*0.2))
            sys_q4.append(mean(sys_qq))

    return [sys_q3,sys_q4]


def funfind(funv,funq3,funq4):
    funw=0.0
    funt=True
    funi=0
    while funt:
        if (funq3[funi]<=funv and funq3[funi+1]>=funv):
            funt=False
            funqa=(funq4[funi]-funq4[funi+1])/(funq3[funi]-funq3[funi+1])
            funqb=funq4[funi]-funqa*funq3[funi]
            funw=funqa*funv+funqb
        else:
            funi+=1
        if (funi==len(funq3)-1):
            funt=False
    return funw


def logL(logL_V,logL_BC):
    return (1.896-0.4*(logL_V+logL_BC))


def logTN(TN_Teff,TN_logL):
    global ZAMS
    # TN_Teff must be in absolute value, not log
    logTZAMS=-9999.9
    ZAMS_T=list(ZAMS[:,0])
    ZAMS_L=list(ZAMS[:,1])
    
    if (TN_logL>=ZAMS_L[0] and TN_logL<=ZAMS_L[-1]):
        TN_i=0
        TN_found=0
        while (TN_found==0):       
            if ((TN_logL>=ZAMS_L[TN_i]) and (TN_logL<=ZAMS_L[TN_i+1])):
                logTZAMS=TN_logL*(ZAMS_T[TN_i+1]-ZAMS_T[TN_i])/(ZAMS_L[TN_i+1]-ZAMS_L[TN_i])+ZAMS_T[TN_i]-(ZAMS_T[TN_i+1]-ZAMS_T[TN_i])/(ZAMS_L[TN_i+1]-ZAMS_L[TN_i])*ZAMS_L[TN_i]
                #logTZAMS=ZAMS_T[TN_i]
                TN_found=1
            elif (TN_i<len(ZAMS_T)-1):
                TN_i+=1
            else:
                TN_found=1
    elif (TN_logL<ZAMS_L[0]):
        logTZAMS=TN_logL*(ZAMS_T[1]-ZAMS_T[0])/(ZAMS_L[1]-ZAMS_L[0])+ZAMS_T[0]-(ZAMS_T[1]-ZAMS_T[0])/(ZAMS_L[1]-ZAMS_L[0])*ZAMS_L[0]
    else:
        logTZAMS=TN_logL*(ZAMS_T[-1]-ZAMS_T[-2])/(ZAMS_L[-1]-ZAMS_L[-2])+ZAMS_T[-2]-(ZAMS_T[-1]-ZAMS_T[-2])/(ZAMS_L[-1]-ZAMS_L[-2])*ZAMS_L[-2]
    
    return log10(TN_Teff)-logTZAMS


def logTZ(TN_Teff,TN_logL):
    global ZAMS
    # TN_Teff must be in absolute value, not log
    logTZAMS=-9999.9
    ZAMS_T=list(ZAMS[:,0])
    ZAMS_L=list(ZAMS[:,1])
    
    if (TN_logL>=ZAMS_L[0] and TN_logL<=ZAMS_L[-1]):
        TN_i=0
        TN_found=0
        while (TN_found==0):       
            if ((TN_logL>=ZAMS_L[TN_i]) and (TN_logL<=ZAMS_L[TN_i+1])):
                logTZAMS=TN_logL*(ZAMS_T[TN_i+1]-ZAMS_T[TN_i])/(ZAMS_L[TN_i+1]-ZAMS_L[TN_i])+ZAMS_T[TN_i]-(ZAMS_T[TN_i+1]-ZAMS_T[TN_i])/(ZAMS_L[TN_i+1]-ZAMS_L[TN_i])*ZAMS_L[TN_i]
                #logTZAMS=ZAMS_T[TN_i]
                TN_found=1
            elif (TN_i<len(ZAMS_T)-1):
                TN_i+=1
            else:
                TN_found=1
    elif (TN_logL<ZAMS_L[0]):
        logTZAMS=TN_logL*(ZAMS_T[1]-ZAMS_T[0])/(ZAMS_L[1]-ZAMS_L[0])+ZAMS_T[0]-(ZAMS_T[1]-ZAMS_T[0])/(ZAMS_L[1]-ZAMS_L[0])*ZAMS_L[0]
    else:
        logTZAMS=TN_logL*(ZAMS_T[-1]-ZAMS_T[-2])/(ZAMS_L[-1]-ZAMS_L[-2])+ZAMS_T[-2]-(ZAMS_T[-1]-ZAMS_T[-2])/(ZAMS_L[-1]-ZAMS_L[-2])*ZAMS_L[-2]
    
    return logTZAMS


def isochrone_grid(grid_age,grid_z,grid_v):
    global isochronesLTN
    global isochronesCMD
    global isochronesZMS
    global isochronesTEF
    global age_values
    global z_values
    if (grid_v=='LTN'):
        grid_iso=isochronesLTN[age_values.index(grid_age)][z_values.index(grid_z)]
        grid_TN=[]
        for grid_i in range(len(array(grid_iso)[:,1])):
            grid_TN.append(logTN(10.0**array(grid_iso)[grid_i,1] , array(grid_iso)[grid_i,0]))
        return [array(grid_iso)[:,0] , array(grid_TN)]
    elif (grid_v=='CMD'):
        grid_iso=isochronesCMD[age_values.index(grid_age)][z_values.index(grid_z)]
        return [array(grid_iso)[:,0] , array(grid_iso)[:,1]]
    elif (grid_v=='ZMS'):
        grid_iso=isochronesZMS[age_values.index(grid_age)][z_values.index(grid_z)]
        grid_TZ=[]
        for grid_i in range(len(array(grid_iso)[:,1])):
            grid_TZ.append(logTZ(10.0**array(grid_iso)[grid_i,1] , array(grid_iso)[grid_i,0]))
        return [array(grid_iso)[:,0] , array(grid_TZ)]
    elif (grid_v=='TEF'):
        grid_iso=isochronesTEF[age_values.index(grid_age)][z_values.index(grid_z)]
        return [array(grid_iso)[:,0] , array(grid_iso)[:,1]]
#-------------------------------------------------------
#-------------------END CALCULATIONS--------------------
#-------------------------------------------------------





#-------------------------------------------------------
#---------------------BINARY CHECK----------------------
#------------------------------------------------------- 
# DO NOT USE !
# DO NOT USE !
# DO NOT USE !
def BinaryJob(dist,clrexc,metal,doname,filt,b,expcor):
    # expcor should be either 0 or 1
    if (expcor==0):
        pass
    else:
        expcor=1

    # transformation from Z to Fe/H
    # metal=metal_transf(0,metal)
    # not necessary anymore

    # transformation factor between E(B-V) and the chosen photometric system colour
    factor_clrexc=clrexc_multiplier(dist,b,filt,expcor)

    # input stellar data
    fdata=loadtxt('clusters/'+str(doname)+'.txt', skiprows=1)
    for i in range(len(fdata)):
        fdata[i][1]=fdata[i][1]-clrexc*factor_clrexc

    lumin0=[]
    tempe0=[]
    ident0=[]
    for i in range(len(fdata)):
        lumintest=( logL( absmag(clrexc,fdata[i][0],dist,filt,b,expcor) , BolCorBV(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) ) )
        tempetest=( logTN( Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) , lumintest ) )
        if (tempetest>-3.0 and tempetest<3.0):
            lumin0.append( lumintest )
            tempe0.append( tempetest )
            ident0.append(i)
    
    tempeB=[]
    luminB=[]
    identB=[]
    for i in range(len(tempe0)):
        if (lumin0[i]<1.0 and lumin0[i]>-1.5):
            tempeB.append(tempe0[i])
            luminB.append(lumin0[i])
            identB.append(ident0[i])

    if (len(luminB)>5):
        binlow=polyfit(luminB,tempeB,3)
        tempeZ=[]
        for i in range(len(tempeB)):
            tempez=0.0
            for j in range(len(binlow)):
                tempez+=binlow[j]*luminB[i]**(len(binlow)-1-j)
            tempeZ.append(tempeB[i]-tempez)
        # skew was supposed to be a measure of the binarity, but it does not work
        #print(skew(tempeZ))
        
        saveh=histogram(tempeZ,bins=10)
        savea=[]
        saveb=[]
        for i in range(len(saveh[0])):
            savea.append(0.5*(saveh[1][i]+saveh[1][i+1]))
            saveb.append(saveh[0][i])
        limrng=savea[saveb.index(max(saveb))]-0.04*0.5

        saveid=[]
        for i in range(len(tempeZ)):
            # if smaller than possibly a binary our an outlier
            if (tempeZ[i]>limrng):
                saveid.append(identB[i])
        for i in range(len(tempe0)):
            if (lumin0[i]>=1.0 or lumin0[i]<=-1.5):
                saveid.append(ident0[i])
        
        return [saveid,0.0]
    else:
        return [[],0.0]
#-------------------------------------------------------
#-----------------END BINARY CHECK----------------------
#------------------------------------------------------- 





#-------------------------------------------------------
#---------------------ISOCHR CHECK----------------------
#------------------------------------------------------- 
def DoJob(dist,clrexc,metal,doname,filt,b,expcor,bincor,binlist):
    # isochrone information from the main body
    global isochronesLTN
    global isochronesCMD
    global isochronesZMS
    global isochronesTEF
    global age_values
    global z_values

    # expcor and bincor should be either 0 or 1 (1 only when starting from ext.maps)
    if (expcor==0):
        pass
    else:
        expcor=1
    if (bincor==0):
        pass
    else:
        bincor=1

    # transformation from Z to Fe/H, based on Pöhnl & Paunzen (2010)
    # metal=metal_transf(0,metal)

    # transformation factor between E(B-V) and the chosen photometric system colour
    factor_clrexc=clrexc_multiplier(dist,b,filt,expcor)

    # input stellar data and transform from CMD space to logL-TN space
    # if bincor=1 then the list of indices will be used for fdata instead of the whole list
    # only lumin0, tempe0 and count0 are necessary for the procedure, the rest is for debugging
    if (bincor==0):
        fdata=loadtxt('clusters/'+str(doname)+'.txt', skiprows=1)
        BpRp=[]
        for i in range(len(fdata)): 
            BpRp.append(0.0+fdata[i][1])
            fdata[i][1]=fdata[i][1]-clrexc*factor_clrexc
        lumin0=[]
        tempe0=[]
        count0=[]
        Gmag0=[]
        Gmag1=[]
        Gmag2=[]
        BpRp0=[]
        BpRp1=[]
        AuxBC=[]
        AuxTT=[]
        AuxTZ=[]
        for i in range(len(fdata)):
            lumintest=( logL( absmag(clrexc,fdata[i][0],dist,filt,b,expcor) , BolCorBV(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) ) )
            tempetest=( logTN( Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) , lumintest ) )
            if (tempetest>-3.0 and tempetest<3.0 and lumintest>-5.0 and lumintest<5.0):
                lumin0.append( lumintest )
                tempe0.append( tempetest )
                count0.append( 1.0 )
                Gmag0.append(absmag(clrexc,fdata[i][0],dist,filt,b,expcor))
                Gmag1.append(fdata[i][0])
                Gmag2.append(absmag(clrexc,fdata[i][0],dist,filt,b,expcor)+5*log10(dist)-5.0)
                BpRp0.append(fdata[i][1])
                BpRp1.append(BpRp[i])
                AuxBC.append(BolCorBV(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt))
                AuxTT.append(Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt))
                AuxTZ.append(logTZ( Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) , lumintest ))
    else:
        fdata=loadtxt(str(doname)+'.txt', skiprows=1)
        BpRp=[]
        for i in binlist: 
            BpRp.append(0.0+fdata[i][1])
            fdata[i][1]=fdata[i][1]-clrexc*factor_clrexc
        lumin0=[]
        tempe0=[]
        count0=[]
        Gmag0=[]
        Gmag1=[]
        Gmag2=[]
        BpRp0=[]
        BpRp1=[]
        AuxBC=[]
        AuxTT=[]
        AuxTZ=[]
        for i in binlist:
            lumintest=( logL( absmag(clrexc,fdata[i][0],dist,filt,b,expcor) , BolCorBV(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) ) )
            tempetest=( logTN( Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) , lumintest ) )
            if (tempetest>-3.0 and tempetest<3.0 and lumintest>-5.0 and lumintest<5.0):
                lumin0.append( lumintest )
                tempe0.append( tempetest )
                count0.append( 1.0 )
                Gmag0.append(absmag(clrexc,fdata[i][0],dist,filt,b,expcor))
                Gmag1.append(fdata[i][0])
                Gmag2.append(absmag(clrexc,fdata[i][0],dist,filt,b,expcor)+5*log10(dist)-5.0)
                BpRp0.append(fdata[i][1])
                BpRp1.append(BpRp[i])
                AuxBC.append(BolCorBV(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt))
                AuxTT.append(Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt))
                AuxTZ.append(logTZ( Teff(fdata[i][1],absmag(clrexc,fdata[i][0],dist,filt,b,expcor),metal,filt) , lumintest ))

    # finding the best fit in isochrone grid
    fitI=-1
    fitJ=-1
    fitX=1.0e16
    for i in range(len(age_values)):
        for j in range(len(z_values)):
            tempe1=[]
            # apply TN systematic corrections, fitting section
            for l in range(len(tempe0)):
                tempe1.append(tempe0[l]+funfind(lumin0[l],temp_corr_grid[i][j][0],temp_corr_grid[i][j][1]))
            fitvalues=LstSqr(lumin0,tempe1,count0,isochronesLTN[i][j][0],isochronesLTN[i][j][1])
            if (fitvalues<fitX):
                fitX=fitvalues
                fitI=0+i
                fitJ=0+j

    # apply TN systematic corrections, results section
    AuxSY=[]
    for l in range(len(tempe0)):
        AuxSY.append(funfind(lumin0[l],temp_corr_grid[fitI][fitJ][0],temp_corr_grid[fitI][fitJ][1]))
        tempe0[l]+=funfind(lumin0[l],temp_corr_grid[fitI][fitJ][0],temp_corr_grid[fitI][fitJ][1])

    return([age_values[fitI],z_values[fitJ],fitX,[tempe0,lumin0,isochronesLTN[fitI][fitJ][1],isochronesLTN[fitI][fitJ][0],AuxSY,isochronesZMS[fitI][fitJ][1],isochronesTEF[fitI][fitJ][1]],[BpRp0,Gmag0,isochronesCMD[fitI][fitJ][1],isochronesCMD[fitI][fitJ][0],BpRp1,Gmag1,Gmag2,AuxBC,AuxTT,AuxTZ]])
#-------------------------------------------------------
#-----------------END ISOCHR CHECK----------------------
#------------------------------------------------------- 




#<<<<---------------------- ------------------- ---------------------->>>>
#<<<<---------------------- ------------------- ---------------------->>>>
#<<<<---------------------- PROGRAM STARTS HERE ---------------------->>>>
#<<<<---------------------- ------------------- ---------------------->>>>
#<<<<---------------------- ------------------- ---------------------->>>>

# only for testing
debugTest=False


# initialisation, loads photometric system, grid spacing and iteration numbers
photosystem,age_step,z_step,Nredd,Niter,redAdj=initialise()


print('Loading list of clusters ...')
try:
    pathfile=os.getcwd()
    os.chdir(pathfile)
    with open("clusters/_complete.txt","r") as f_data:
        all_data=[x.split() for x in f_data.readlines()]
        dataini=array([list(map(str,x)) for x in all_data[1:]])
except FileNotFoundError:
    pathfile=os.path.dirname(__file__)
    os.chdir(pathfile)
    with open("clusters/_complete.txt","r") as f_data:
        all_data=[x.split() for x in f_data.readlines()]
        dataini=array([list(map(str,x)) for x in all_data[1:]])
clust_list=[]
ini_D_list=[]
ini_E_list=[]
exp_A_list=[]
exp_Z_list=[]
par_b_list=[]
for i in range(len(dataini)):
    clust_list.append(dataini[i][0])
    ini_D_list.append(float(dataini[i][3]))
    ini_E_list.append(float(dataini[i][4]))
    par_b_list.append(float(dataini[i][1]))
    exp_A_list.append(-999)
    exp_Z_list.append(-999)

print('Preparing isochrone grid ...')

# load ZAMS for the calculations
with open("ZAMS_014.txt","r") as f_data:
    all_data=[x.split() for x in f_data.readlines()]
    ZAMS=array([list(map(float,x)) for x in all_data[0:]])

# preparing isochrone grid
with open("isochrones"+photosystem+".txt","r") as f_data:
    all_data=[x.split() for x in f_data.readlines()]
    isochrones_complete=array([list(map(str,x)) for x in all_data[13:]])
f_data=[]
all_data=[]

# first, create a list of ages and metallicities available
age_values=[6.6]
# age_step=0.2
age_last=10.0
while (age_values[-1]!=age_last):
    age_values.append(round(age_values[-1]+age_step , 1))
z_values=[0.005]
# z_step=0.005
z_last=0.040
while (z_values[-1]!=z_last):
    z_values.append(round(z_values[-1]+z_step , 3))

# create the grid, using age for rows and metallicity for columns
isochronesLTN=[]
isochronesCMD=[]
isochronesZMS=[]
isochronesTEF=[]
for i in range(len(age_values)):
    isohelp1=[]
    isohelp2=[]
    isohelp3=[]
    isohelp4=[]
    for j in range(len(z_values)):
        isohelp1.append([])
        isohelp2.append([])
        isohelp3.append([])
        isohelp4.append([])
    isochronesLTN.append(isohelp1)
    isochronesCMD.append(isohelp2)
    isochronesZMS.append(isohelp3)
    isochronesTEF.append(isohelp4)
    
# fill in the grid from the loaded data (isochrones_complete)
for i in range(len(isochrones_complete)):
    if (age_values.count(round(float(isochrones_complete[i][2]),1))==1 and z_values.count(round(float(isochrones_complete[i][0]),3))==1):
        age_idx=age_values.index(round(float(isochrones_complete[i][2]),1))
        z_idx=z_values.index(round(float(isochrones_complete[i][0]),3))
        isochronesLTN[age_idx][z_idx].append([float(isochrones_complete[i][6]) , float(isochrones_complete[i][7])])
        isochronesZMS[age_idx][z_idx].append([float(isochrones_complete[i][6]) , float(isochrones_complete[i][7])])
        isochronesTEF[age_idx][z_idx].append([float(isochrones_complete[i][6]) , float(isochrones_complete[i][7])])
        if (photosystem=='G'):
            isochronesCMD[age_idx][z_idx].append([float(isochrones_complete[i][28]) , (float(isochrones_complete[i][29])-float(isochrones_complete[i][30]))])
        elif (photosystem=='2'):
            isochronesCMD[age_idx][z_idx].append([float(isochrones_complete[i][28]) , (float(isochrones_complete[i][28])-float(isochrones_complete[i][30]))])
        elif (photosystem=='J'):
            isochronesCMD[age_idx][z_idx].append([float(isochrones_complete[i][30]) , (float(isochrones_complete[i][29])-float(isochrones_complete[i][30]))])
isochrones_complete=[]

# transform isochrones to the normalised grid
for main_i in range(len(age_values)):
    for main_j in range(len(z_values)):
        isochronesLTN[main_i][main_j]=isochrone_grid(age_values[main_i],z_values[main_j],'LTN')
        isochronesCMD[main_i][main_j]=isochrone_grid(age_values[main_i],z_values[main_j],'CMD')
        isochronesZMS[main_i][main_j]=isochrone_grid(age_values[main_i],z_values[main_j],'ZMS')
        isochronesTEF[main_i][main_j]=isochrone_grid(age_values[main_i],z_values[main_j],'TEF')

# prepare grid for corrected systematics for temperature calibration
temp_corr_grid=[]
for main_i in range(len(age_values)):
    temp_corr_help=[]
    for main_j in range(len(z_values)):
        temp_corr_help.append(sys_temp(age_values[main_i],z_values[main_j],photosystem))
    temp_corr_grid.append(temp_corr_help)

print('Starting main procedure ...\n')
start_time = time.time()
for cc in range(len(clust_list)):
    try:
        start_time_1 = time.time()
        # starting values (colour excess is E(B-V), this is recalculated to whatever system
        # using the factor_clrexc)
        iniD=ini_D_list[cc]
        iniE=ini_E_list[cc]
        iniZ=0.014
        clustname=clust_list[cc]+'_'+photosystem

        # numbX = number of iterations of given quantity
        # rangX = range of itereations centred around the ini value
        # valX = current value in iteration
        numbD=1
        rangD=0.8*iniD
        if (iniE<0.0 or Nredd==0):
            numbE=10
            rangE=[0.010,0.040,0.080,0.125,0.250,0.500,0.750,1.000,1.500,2.000]
        elif (Nredd==1):
            numbE=0+Nredd
            rangE=[iniE]
        else:
            numbE=0+Nredd
            rangE=[]
            for main_i in range(numbE):
                rangE.append(iniE-redAdj*iniE+2*redAdj*iniE*main_i/float(numbE-1))

        # start calculating fits for the grid of parameters, assume no binaries
        valD=iniD-0.0*rangD
        result_grid=[]
        for main_i in range(numbD):
            for main_j in range(numbE):
                valE=rangE[main_j]
                final_z=0.000+iniZ
                res_ini=0.000
                check_limit=0
                while ((round(res_ini,3)!=round(final_z,3)) and check_limit<Niter):
                    res_ini=0.000+final_z
                    final_age,final_z,final_fit,final_data,final_data2 = DoJob(valD,valE,res_ini,clustname,photosystem,par_b_list[cc],0,0,[])
                    check_limit+=1
                result_grid.append([valD,valE,final_age,final_z,final_fit,final_data,final_data2,check_limit])
                print(check_limit)
            valD+=rangD/float(numbD-1+1)

        # results for all reddening values are sorted are written in a logfile
        print(clust_list[cc])
        result_grid=array(result_grid)
        sorted_result_grid=result_grid[result_grid[:,4].argsort()]
        print('\n')
        fsave = open('finished/_logfile.txt', "a+")
        fsave.write('Cluster: %s \n' % (clust_list[cc]))
        fsave.write('Inputs: %s ; %f ; %f ; %d ; %f ; %d \n' % (photosystem,age_step,z_step,Nredd,redAdj,Niter))
        fsave.write('------------------------------------------------------------------------ \n')
        for main_i in range(len(sorted_result_grid)):
            fsave.write('Parameters:   %d   %.3f   %.1f   %.3f  .....  fit:%.7f   iter:%d \n' % (int(sorted_result_grid[main_i][0]),sorted_result_grid[main_i][1],sorted_result_grid[main_i][2],sorted_result_grid[main_i][3],sorted_result_grid[main_i][4],sorted_result_grid[main_i][7]))
        fsave.write('# FINISHED (%.3f min) \n\n' % ((time.time() - start_time_1)/60.0))
        fsave.close()
        
        # the three best results are plotted in CMD and LTN diagrams
        # if debug mode is on, additional data for individual points are returned for the three fits
        fitcurves=['r-','g--','b:']
        fitpoints=['ko','ko','ko']
        plt.figure(figsize=(12,6))

        for main_i in range(len(fitcurves)):
            try:
                plt.subplot(2,3,main_i+1)
                plt.plot(sorted_result_grid[main_i][5][0],sorted_result_grid[main_i][5][1],fitpoints[main_i],ms=4.0,alpha=0.2)
                plt.plot(sorted_result_grid[main_i][5][2],sorted_result_grid[main_i][5][3],fitcurves[main_i],alpha=0.6)
                if (main_i==0):
                    plt.xlabel('TN')
                    plt.ylabel('log L')
                else:
                    plt.xlabel('TN')
                #plt.xlim(-0.8,0.2)
                #plt.ylim(-2.5,4.0)
                plt.xlim(  min(sorted_result_grid[main_i][5][0]) - 0.25*(max(sorted_result_grid[main_i][5][0])-min(sorted_result_grid[main_i][5][0]))  ,  max(sorted_result_grid[main_i][5][0]) + 0.25*(max(sorted_result_grid[main_i][5][0])-min(sorted_result_grid[main_i][5][0]))  )
                plt.ylim(  min(sorted_result_grid[main_i][5][1]) - 0.05*(max(sorted_result_grid[main_i][5][1])-min(sorted_result_grid[main_i][5][1]))  ,  max(sorted_result_grid[main_i][5][1]) + 0.05*(max(sorted_result_grid[main_i][5][1])-min(sorted_result_grid[main_i][5][1]))  )
                plt.title('%d ; %.3f ; %.1f ; %.3f ... %.6f' % (int(sorted_result_grid[main_i][0]),sorted_result_grid[main_i][1],sorted_result_grid[main_i][2],sorted_result_grid[main_i][3],sorted_result_grid[main_i][4]))
                plt.locator_params(axis='x',nbins=7)

                if (debugTest):
                    fsave = open('finished/'+str(clustname)+'_aux_R'+str(main_i+1)+'_isoch.txt', "w+")
                    fsave.write('Cluster: %s \n' % (clustname))
                    fsave.write('Parameters: %d ; %.3f ; %.1f ; %.3f ... %.6f \n' % (int(sorted_result_grid[main_i][0]),sorted_result_grid[main_i][1],sorted_result_grid[main_i][2],sorted_result_grid[main_i][3],sorted_result_grid[main_i][4]))
                    fsave.write('Inputs: %s ; %f ; %f ; %d ; %f ; %d \n' % (photosystem,age_step,z_step,Nredd,redAdj,Niter))
                    fsave.write('---------------------------------------------------------------- \n')
                    fsave.write('%7s %7s %7s %7s \n' % ('T_eff','T_ZAMS','T_N','logL'))
                    for main_j in range(len(sorted_result_grid[main_i][5][2])):
                        fsave.write('%7.4f %7.4f %7.4f %7.3f \n' % (sorted_result_grid[main_i][5][6][main_j],sorted_result_grid[main_i][5][5][main_j],sorted_result_grid[main_i][5][2][main_j],sorted_result_grid[main_i][5][3][main_j]))
                    fsave.close()
            except IndexError:
                if ((Nredd==1 and main_i>0)==False):
                    print('err')

            try:
                plt.subplot(2,3,main_i+4)
                plt.plot(sorted_result_grid[main_i][6][0],sorted_result_grid[main_i][6][1],fitpoints[main_i],ms=4.0,alpha=0.2)
                plt.plot(sorted_result_grid[main_i][6][2],sorted_result_grid[main_i][6][3],fitcurves[main_i],alpha=0.6)
                if (main_i==0):
                    if (photosystem=='G'):
                        plt.xlabel('(BP-RP)_0 [mag]')
                        plt.ylabel('M_G [mag]')
                    elif (photosystem=='2'):
                        plt.xlabel('(J-Ks)_0 [mag]')
                        plt.ylabel('M_J [mag]')
                    elif (photosystem=='J'):
                        plt.xlabel('(B-V)_0 [mag]')
                        plt.ylabel('M_V [mag]')
                else:
                    if (photosystem=='G'):
                        plt.xlabel('(BP-RP)_0 [mag]')
                    elif (photosystem=='2'):
                        plt.xlabel('(J-Ks)_0 [mag]')
                    elif (photosystem=='J'):
                        plt.xlabel('(B-V)_0 [mag]')
                # plt.xlim(-0.8,4.0)
                # plt.ylim(15,-5)
                plt.xlim(  min(sorted_result_grid[main_i][6][0]) - 0.10*(max(sorted_result_grid[main_i][6][0])-min(sorted_result_grid[main_i][6][0]))  ,  max(sorted_result_grid[main_i][6][0]) + 0.10*(max(sorted_result_grid[main_i][6][0])-min(sorted_result_grid[main_i][6][0]))  )
                plt.ylim(  max(sorted_result_grid[main_i][6][1]) + 0.10*(max(sorted_result_grid[main_i][6][1])-min(sorted_result_grid[main_i][6][1]))  ,  min(sorted_result_grid[main_i][6][1]) - 0.10*(max(sorted_result_grid[main_i][6][1])-min(sorted_result_grid[main_i][6][1]))  )
                #plt.gca().invert_yaxis()
                #plt.title('%d ; %.3f ; %.1f ; %.3f ... %.6f' % (int(sorted_result_grid[main_i][0]),sorted_result_grid[main_i][1],sorted_result_grid[main_i][2],sorted_result_grid[main_i][3],sorted_result_grid[main_i][4]))

                if (debugTest):
                    fsave = open('finished/'+str(clustname)+'_aux_R'+str(main_i+1)+'_clust.txt', "w+")
                    fsave.write('Cluster: %s \n' % (clustname))
                    fsave.write('Parameters: %d ; %.3f ; %.1f ; %.3f ... %.6f \n' % (int(sorted_result_grid[main_i][0]),sorted_result_grid[main_i][1],sorted_result_grid[main_i][2],sorted_result_grid[main_i][3],sorted_result_grid[main_i][4]))
                    fsave.write('Inputs: %s ; %f ; %f ; %d ; %f ; %d \n' % (photosystem,age_step,z_step,Nredd,redAdj,Niter))
                    fsave.write('---------------------------------------------------------------- \n')
                    fsave.write('%6s %6s %6s %6s %6s %6s %7s %7s %7s %7s %7s \n' % ('Color','Color0','Mag','Mag0','AbsMag','BC_Mag','T_eff','T_ZAMS','T_corr','T_N','logL'))
                    for main_j in range(len(sorted_result_grid[main_i][5][0])):
                        fsave.write('%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.4f %7.4f %7.4f %7.4f %7.3f \n' % (sorted_result_grid[main_i][6][4][main_j],sorted_result_grid[main_i][6][0][main_j],sorted_result_grid[main_i][6][5][main_j],sorted_result_grid[main_i][6][6][main_j],sorted_result_grid[main_i][6][1][main_j],sorted_result_grid[main_i][6][7][main_j],log10(sorted_result_grid[main_i][6][8][main_j]),sorted_result_grid[main_i][6][9][main_j],sorted_result_grid[main_i][5][4][main_j],sorted_result_grid[main_i][5][0][main_j],sorted_result_grid[main_i][5][1][main_j]))
                    fsave.close()
            except IndexError:
                if ((Nredd==1 and main_i>0)==False):
                    print('err')

        plt.tight_layout()
        plt.savefig('finished/'+str(clustname)+'.png',dpi=300,bbox_inches="tight")
        plt.close()

    except OSError:
        # exception can be encountered if the list of clusters does not match the provided data files
        # names of the clusters should coincide with the names of the data files (up to the photo.system designation)
        print('no file: %s' % (clustname))

isochrones=[]
print('\n')
print('Finished!')
print("Runtime: %s min" % ((time.time() - start_time)/60.0))