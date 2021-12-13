def metal_transf(transf_way,transf_numb):
    transf_feh=[-0.729,-0.525,-0.387,-0.282,-0.224,-0.149,-0.086,-0.030,0.019,0.077,0.116,0.152,0.185,0.225,0.253,0.288,0.312,0.343,0.371]
    transf_z=[0.004,0.006,0.008,0.010,0.012,0.014,0.016,0.018,0.020,0.022,0.024,0.026,0.028,0.030,0.032,0.034,0.036,0.038,0.040]
    transf_fin=0.0
    if (transf_way==0): #z to Fe/H
        transf_found=0
        transf_j=0
        while (transf_found==0):
            if (transf_numb<=transf_z[transf_j]):
                transf_fin=transf_feh[transf_j]
                transf_found=1
            else:
                transf_j+=1
    elif (transf_way==1): #Fe/H to z
        transf_found=0
        transf_j=0
        while (transf_found==0):
            if (transf_numb<=transf_feh[transf_j]):
                transf_fin=transf_z[transf_j]
                transf_found=1
            else:
                transf_j+=1
    return transf_fin