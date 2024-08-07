# -*- coding: utf-8 -*-
#Spyder Editor
#
#Python version of Caroline's CS model'
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
   
# initialize model setting
   
MAXPROF=26      # MAXIMUM NUMBER OF PROFILES
MAXNOUR=10      # MAXIMUM NUMBER OF NOURISHMENT
MAXTIME=65000   # MAXIMUM NUMBER OF TIME STEPS 

# start with subroutines
# %%

def FALVEL(TEMP,D50):
    print('computing fall velocity')
    SRATIO=2.65
    NYVAL = np.zeros(9)
    TETA = np.zeros(9)
    NYVAL[0]=1.79
    NYVAL[1]=1.51
    NYVAL[2]=1.31
    NYVAL[3]=1.14
    NYVAL[4]=1.00
    NYVAL[5]=0.894
    NYVAL[6]=0.799
    NYVAL[7]=0.724
    NYVAL[8]=0.658
    NYVAL=NYVAL*1E-6
    
    for I in range(9):
        TETA[I]=float(I)*5.0
 
        
    if (TEMP<TETA[0]) or (TEMP>TETA[8]):
      print('TEMPERATURE FOR VISCOSITY CALCULATION IS OUT OF RANGE')
      print(TEMP)
      
# linear interpolation
    NY = np.interp(TEMP,TETA,NYVAL)     
    B=(SRATIO-1)*9.81*D50**3/NY**2
    
    if (B<39):
        WFALL=(SRATIO-1)*9.81*D50**2/18/NY
    elif (B>39 and B<1E4):
        WFALL=((SRATIO-1)*9.81)**0.7*D50**1.1/6/NY**0.4
    else:    
        WFALL=((SRATIO-1)*9.81*D50/0.91)**0.5
    return WFALL
    
# %%   
import math
def subar(DT, YS, YL, S, HO, T, DFOOT, D50, R, VBEACH, CF, WSPEED, WDIR, KW, CWRATE, CIMPACT, QS, SWL,
          BETAF, YG, QWNET, QSLR, QLS, VNOURB, QNOUR, VW, QL, K, A, APRIM, ASFN):

    G = 9.81
    PI = 3.14159

    X = 2 * VBEACH / DFOOT * (1 - SWL / DFOOT)
    LO = 1.5613 * (0.95 * T) ** 2
    R = 0.158 * math.sqrt(HO * LO)
    RPRIM = R * math.exp(-2. * CF * X) + (DFOOT - SWL) * (1. - math.exp(-2. * CF * X))

    if RPRIM + SWL > DFOOT:
        R = RPRIM

    VWNEW = VW + (QS - QWNET - A * (QSLR - QLS) + ASFN * QNOUR) * DT + APRIM * VNOURB

   
    RHOS = 2650.
    RHOA = 1.2
    D50REF = 0.25e-3
    AWIND = 0.1
    P = 0.4
    Z = 10
    KAPPA = 0.4

    WCRIT = AWIND * math.sqrt((RHOS - RHOA) / RHOA * G * D50)

    Z0 = 0.081 * math.log10(D50 * 10 ** 3 / 0.18)
 

    WSHEAR = WSPEED * KAPPA * 1 / math.log(Z / Z0)

    if WSHEAR > WCRIT and WSPEED > 7:
        MWE = KW * math.sqrt(D50 / D50REF) * RHOA * WSHEAR ** 2. / G * (WSHEAR - WCRIT)
    else:
        MWE = 0

    QWE = MWE / RHOS * (1 - P)

    YR = YS + (1 - (R + SWL) / DFOOT) * (YG - YS)

    BDRY = YR - YS

    if R + SWL >= DFOOT or YR < YS:
        BDRY = 0

    F = BDRY / abs(math.cos(WDIR))
    QW = QWE * (1. - math.exp(-CWRATE * F))
    QWNET = QW * math.cos(WDIR)



    if QWNET * DT > VW:
        QWNET = VW / DT

# introduce dynamic reduction based on Duran and Moore [2013] Ad
# BETA IS SET AT 4.2 HERE
    QWNET = QWNET*max(1.0-4.2*S/(YG-YS),0.0)    
#    QWNET = QWNET*max(1.0-4.2*S/(50),0.0)   
#    QWNET = 0.0    

    if VWNEW <= 0:
        VWNEW = 0
        QWNET = 0
        QWE = 0
        F = 0
        BDRY = 0
        WSHEAR = 0
        WCRIT = 0
        QW = 0
         
    VW = VWNEW

    AOVER = 3.0
    ALFA = 0.0

    if R + SWL > DFOOT:
        if R + SWL <= DFOOT + S:
            QD = 4. * CIMPACT * (R + SWL - DFOOT) ** 2 / T
        else:
            print('OVERWASH')
            QD = 4. * CIMPACT * (R + SWL - DFOOT) * S / T
            ALFA = ((R + SWL - DFOOT) / S - 1.) / AOVER
    else:
        QD = 0.0

    QS = QD / (1 + ALFA)
    QL = QD * ALFA / (1 + ALFA)
    
    
    return QWNET,QS,QL,VW,R    

# TEST CALL
# DT = 10800
# YS = 147.5
# YL = 99.8
# S = 4.16
# HO = 2.16
# T = 11.33
# DFOOT = 3.0
# D50 = 0.0002
# R = 0.0
# VBEACH = 446
# CF = 0.03
# WSPEED = 14.4
# WDIR = 1.0
# KW = 9.5
# CWRATE = 0.1
# CIMPACT = 0.0006
# QS = 8.6E-06
# SWL=3.01
# BETAF = 0.1
# YG = 387
# QWNET =0.
# QSLR= 4.1E-08
# QLS= 0.0
# VNOURB = 0.0
# QNOUR = 0.0
# VW = 0.0
# QL = 0.0
# K = 0
# A = 0.6
# APRIM = 0.4
# ASFN = 0.2

# QWNET,QS,QL,VW,R=  subar(DT, YS, YL, S, HO, T, DFOOT, 
#                                          D50, R, VBEACH, CF, WSPEED, WDIR, 
#                                          KW, CWRATE, CIMPACT, QS, SWL,
#                                          BETAF, YG, QWNET, QSLR, QLS, 
#                                          VNOURB, QNOUR, VW, QL, K, A, APRIM, ASFN)


#%%
        
def subaq(DT, VOLBAR, HO, T, WFALL, SLR, BARNOUR, CB, BACTIVE):
    # Empirical coefficients
    ACOEFF = 1.333
    BCOEFF = 1.0
    DCOEFF = 0.56e-6
    MCOEFF = -0.5

    # Dean parameter
    DEAN = HO / WFALL / T

    # Deepwater wave steepness
    LO = 1.56 * T**2
    STP = HO / LO

    # Equilibrium bar volume
    VOLBAREQ = CB * DEAN**ACOEFF * STP**BCOEFF * LO**2

    # Response coefficient for bar volume
    LAMBDA2 = DCOEFF * DEAN**MCOEFF

    if VOLBAREQ - VOLBAR > 0:
        LAMBDA = LAMBDA2
    else:
        LAMBDA = 0.30 * LAMBDA2  # You can adjust this value as needed

    # Calculate bar volume change
    DVOL = (VOLBAREQ - VOLBAR) * (1. - math.exp(-LAMBDA * DT))

    # Transport to (from) bar
    QBAR = DVOL / DT

    # Transport onshore from shoreface nourishment
    NOURCOEFF = 0.1
    DVOLNOUR = (BARNOUR) * (1. - math.exp(-NOURCOEFF * LAMBDA * DT))
    QNOUR = DVOLNOUR / DT

    # Bruun rule transport due to sea level rise
    QSLR = SLR * BACTIVE

    return QBAR, QNOUR, QSLR, DVOL

# # Example usage:
# DT = 1.0  # Replace with the actual value
# VOLBAR = 20.0  # Replace with the actual value
# HO = 1.0  # Replace with the actual value
# T = 5.0  # Replace with the actual value
# WFALL = 0.02  # Replace with the actual value
# CB = 1.0  # Replace with the actual value
# BARNOUR = 1.0  # Replace with the actual value
# SLR = 0.0  # Replace with the actual value

# QBAR, QNOUR, QSLR, DVOL = subaq(DT, VOLBAR, HO, T, WFALL, SLR, BARNOUR, CB)

# # Use the calculated values as needed
# print(f"QBAR: {QBAR}, QNOUR: {QNOUR}, QSLR: {QSLR}, DVOL: {DVOL}")

# %%

def nour(i, iv, volat, volad, volab, volrb, maxnour):
    vnourd, vnourb, vnourbar = 0.0, 0.0, 0.0
    if i == int(volat[int(iv) - 1]):
        vnourd = volad[int(iv) - 1]
        vnourb = volab[int(iv) - 1]
        vnourbar = volrb[int(iv) - 1]
        print('beach nourishment ', {vnourb})
        print('bar nourishment ', {vnourbar})
        if iv<maxnour:    
            iv = iv+ 1

    return vnourd, vnourb, vnourbar, iv

# # Example usage:
#i = 1001  # Replace with the actual value
#iv = 2  # Replace with the actual value
#maxnour = 2  # Replace with the actual value

# volat = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Replace with the actual values
# volad = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # Replace with the actual values
# volab = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]  # Replace with the actual values
# volrb = [2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]  # Replace with the actual values

# vnourd, vnourb, vnourbar, iv = nour(i, iv, volat, volad, volab, volrb, maxnour)

# # Use the calculated values as needed
# print(f"VNOURD: {vnourd}, VNOURB: {vnourb}, VNOURBAR: {vnourbar}, IV: {iv}")

# %%

def bar(volbar, qbar, dt, vnourbar, barnour, qnour):
    volbar_new = volbar + qbar * dt
    barnour_new = barnour + vnourbar - qnour * dt

    volbar = volbar_new
    barnour = barnour_new

    if volbar < 0:
        print('BAR VOLUME EXCEEDED')

    return volbar, barnour

# # Example usage:
# volbar = 2.0  # Replace with the actual value
# qbar = 0.1  # Replace with the actual value
# dt = 1.0  # Replace with the actual value
# vnourbar = 0.05  # Replace with the actual value
# barnour = 1.0  # Replace with the actual value
# qnour = 0.02  # Replace with the actual value

# volbar, barnour = bar(volbar, qbar, dt, vnourbar, barnour, qnour)

# # Use the updated values as needed
# print(f"Updated VOLBAR: {volbar}, Updated BARNOUR: {barnour}")

#  

# %%
def duneeros(vramp, ysprim, betas, betal, s, erod, yl, ylprim, dfoot):
#    print('in duneeros')
    deltavramp = erod
    vrampnew = vramp + deltavramp

    if vrampnew > 0:
        srampnew = (2 * betas * vrampnew)**0.5
        ysnew = srampnew / betas + ysprim
        ysprimnew = ysprim
        ylprimnew = ylprim
        ylnew = yl
        snew = s

    elif vrampnew < 0:
        ysprimnew = ysprim + vrampnew / s

        if ysprimnew < ylprim:
            wtransprest = (ysprimnew - ylprim) * s
            ylprimnew = (1 / betal * ((ylprim - yl) * s + 2 * wtransprest))**0.5 + yl
            snew = (ylprimnew - yl) * betal
            ysprimnew = ylprimnew
        else:
            ylprimnew = ylprim
            snew = s

        ysnew = ysprimnew
        ylnew = yl
        vrampnew = 0.0
        srampnew = 0.0


    elif vrampnew == 0:
        ysnew = ysprim
        ysprimnew = ysprim
        ylprimnew = ylprim
        ylnew = yl
        srampnew = 0.0
        snew = s


    return vrampnew, ysnew, ysprimnew, snew, ylnew, ylprimnew, srampnew
#
#test values
# vramp = 18.11
# ysprim = 138.7
# betas = 0.47
# betal = 0.26
# s = 4.16
# erod = -24.0
# yl = 99.88
# ylprim = 115.9
# #
# vrampnew, ysnew, ysprimnew, snew, ylnew, ylprimnew, srampnew = duneeros(vramp, ysprim, betas, betal, s, erod, yl, ylprim)

 # %%

def beach(qwnet, qs, qbar, dclos, db, slr, bactive, dt,
          dfoot, vbeach, vbeachtot, betaf, yl, ys,
          vnourb, qls, qslr, vramp, ysprim, betas,
          s, ylprim, betal, al, bl, qnour, sramp):
    """
    Beach subroutine in Python
    """
    # Constants
    pi = 3.14159

    # Calculate new total beach volume
    vbeachtotnew = vbeachtot + (-qwnet + qs - qbar + qls - qslr + qnour) * dt + vnourb
 
    # Check if vbeachtotnew is smaller then vbeachtotmin -> DUNE EROSION
    vbeachtotmin = ys * (dclos + dfoot) + 0.5 * dfoot**2 / betaf + dfoot / betaf * dclos

    if vbeachtotnew < vbeachtotmin:
        erod = vbeachtotnew - vbeachtotmin
#        vramp_new, ys_new, ysprim_new, s_new, yl_new, ylprim_new, sramp_new = duneeros(
#            vramp, ysprim, betas, erod, s, ys, ysprim, s, yl, ylprim, betal)
        vramp_new, ys_new, ysprim_new, s_new, yl_new, ylprim_new, sramp_new = duneeros(vramp, ysprim, betas, betal, s, erod, yl, ylprim)


        # Update parameters
        vramp = vramp_new
        ys = ys_new
        yl = yl_new
        ysprim = ysprim_new
        ylprim = ylprim_new
        s = s_new
        sramp = sramp_new

        vbeachtotnew = vbeachtotmin

    # Calculate new shoreline coordinate
    yg_new = ys + (al * (vbeachtotnew - ys * (dfoot + dclos)) + bl) / (1 + al * dclos)

    # Calculate db and vbeachnew
    db = (vbeachtotnew - ys * (dclos + dfoot) - (yg_new - ys) * dclos) / (yg_new - ys)
    vbeachnew = (yg_new - ys) * db

    # Update vbeach parameters
    vbeachtot = vbeachtotnew
    vbeach = vbeachnew
    yg = yg_new

    return vramp, ys, ysprim, s, yl, ylprim, sramp, vbeachtot,vbeach,yg
#
# vramp, ys, ysprim, s, yl, ylprim, sramp, vbeachtot,vbeach,yg=beach(qwnet, qs, qbar, dclos, db, slr, bactive, dt,
#          dfoot, vbeach, vbeachtot, betaf, yl, ys,
#          vnourb, qls, qslr, vramp, ysprim, betas,
#          s, ylprim, betal, al, bl, qnour)
#

#%%
def sedbud2_1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp):    
# dune evolution for positive sediment budget in case of triangular shape
#                    print('in sebud2_1')
                    ylprim = ysprim
                    if s < smax:
                        snew = math.sqrt(((ys - yl) * s + 2 * wtransp) / (1 / betas + 1 / betal))
                        if snew > smax:
                            s = smax
                            wtransp = 0.5 * (1 / betal + 1 / betas) * (snew ** 2 - s ** 2)
                            snew = s
                            srampnew = snew
                            ylnew = yl
                            ylprimnew = ylprim
                            ysprimnew = ylprimnew + wtransp / snew
                            ysnew = ysprimnew + snew / betas
                            
                        ylprimnew = ylprim   # dune increases in height
                        ysprimnew = ysprim
                        ysnew = ysprimnew + snew / betas
                        ylnew = ylprimnew - snew / betal
                        srampnew = snew
                        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                    
                    elif s >= smax:  # dune builds out to sea
                        snew = s
                        srampnew = snew
                        ylnew = yl
                        ylprimnew = ylprim
                        ysprimnew = ylprimnew + wtransp / snew
                        ysnew = ysprimnew + snew / betas
                        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
    
                    return snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew
                
# %%
def sedbud2_2(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp):  
# dune evolution for positive sediment budget in case of trapezoidal shape
#                    print('in sebud2_2')
                    if s < smax:
                        vtrimax = 0.5 * (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas) # max triangular volume above s

                        if (s + 2 * vtrimax / (ysprim - ylprim)) <= smax:
                            vmax = vtrimax
                        else:
                            vmax = (smax - s) * (ysprim - ylprim - 0.5 * (smax - s) * (1 / betal + 1 / betas))

                        if wtransp <= vmax:
                            snew = s + 2 * vtrimax / (ysprim - ylprim) - np.sqrt(2 * (vtrimax - wtransp) / (1 / betal + 1 / betas))
                            srampnew = snew
                            ylnew = yl
                            ylprimnew = ylprim + (snew - s) / betal
                            ysprimnew = ysprim - (snew - s) / betas
                            ysnew = ys
                            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                        elif wtransp > vmax:
                            wrest = wtransp - vmax
                            wtransp = wrest
                            if vmax < vtrimax:
                                snew = smax
                                ylprimnew = ylprim + (snew - s) / betal
                                ysprimnew = ysprim - (snew - s) / betas
                                ylnew = yl
                                ysprim = ysprimnew
                                ysprimnew = ysprim + wtransp / snew
                                ysnew = ysprimnew + snew / betas
                                vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                                srampnew = (ysnew - ysprimnew)*betas    # added by Ad
                            else:
                                snew = s + (ysprim - ylprim) / (1 / betas + 1 / betal)
                                ylprimnew = ylprim + (snew - s) / betal
                                ysprimnew = ylprimnew
                                s = snew
                                ylprim = ylprimnew
                                ysprim = ysprimnew
                                ylnew = yl  # added by Ad
                                ysnew = ys  # added by Ad
                                vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas # added by Ad
                                srampnew = (ysnew - ysprimnew)*betas    # added by Ad
                    elif s >= smax:   # seaward moving dune front
                        snew = s
                        ylnew = yl
                        ylprimnew = ylprim
                        ysprimnew = ysprim + wtransp / snew
                        ysnew = ys+ wtransp / snew
                        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas  
                        srampnew = (ysnew - ysprimnew)*betas   # added by Ad                              
#
                    return snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew,wtransp    
#
#
# vramp = 19.67
# ys=149.228
# ysprim = 140.08
# betas = 0.47  
# s = 4.3
# yl = 100
# ylprim = 115.54
# betal = 0.26         
# smax = 4.3
# wtransp = 0.042814434042554694

# snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew = sedbud2_2(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)

     

# %%
def sedbud1_0(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp):    
# dune evolution for stable sediment budget
               print('in sebud1, wtransp', {wtransp})
               if ysprim > ylprim:
                   vtrimax = 0.5 * (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas)
                   if s + 2 * vtrimax / (ysprim - ylprim) <= smax:
                       vmax = vtrimax
                   elif s + 2 * vtrimax / (ysprim - ylprim) > smax:    
                       vmax = (smax - s) * (ysprim - ylprim - 0.5 * (smax - s) * (1 / betal + 1 / betas))
#                   else:  # original code
#                       vmax = (smax - s) * (ysprim - ylprim - 0.5 * (smax - s) * (1 / betal + 1 / betas))
                   if wtransp <= vmax:
                       snew = s + 2 * vtrimax / (ysprim - ylprim) - math.sqrt(2 * (vtrimax - wtransp) / (1 / betal + 1 / betas))
                       srampnew = snew
                       ylnew = yl
                       ylprimnew = ylprim + (snew - s) / betal
                       ysprimnew = ysprim - (snew - s) / betas
                       ysnew = ys
                       vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                   elif wtransp > vmax:
                       wrest = wtransp - vmax
                       wtransp = wrest
                       if vmax < vtrimax:
                           snew = smax
                           ylprimnew = ylprim + (snew - s) / betal
                           ysprimnew = ysprim - (snew - s) / betas
                           ylnew = yl
                           ysprim = ysprimnew
                           ysprimnew = ysprim + wtransp / snew
                           ysnew = ysprimnew + snew / betas
                           vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                           srampnew = snew  # added by Ad
                       else:
                           snew = s + (ysprim - ylprim) / (1 / betas + 1 / betal)
                           ylprimnew = ylprim + (snew - s) / betal
                           ysprimnew = ylprimnew
                           s = snew
                           ylprim = ylprimnew
                           ysprim = ysprimnew
                           srampnew = snew  # added by Ad
                           
               return snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew
           
#
# snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew=sedbud1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)
#            
# %%

def sedbud1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp):    
# dune evolution for stable sediment budget
#    print('in sebud1, wtransp', {wtransp})
    if ysprim > ylprim:  # TRAPEZOIDAL SHAPE, DISTRIBUTE SEDIMENT 50/50 ON CREST AND SEAWARD SLOPE
        # CHECK AMOUNT OF SEDIMENT REQUIRED TO ATTAIN TRIANGULAR SHAPE
        vtrimax = 0.5 * (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas)

        if wtransp / 2 <= vtrimax:
            snew = s + (ysprim - ylprim) / (1 / betal + 1 / betas) - (
                    (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas) ** 2 - wtransp / (1 / betas + 1 / betal)) ** 0.5
            ylprimnew = ylprim + (snew - s) / betal
            ysprimnew = ysprim - (snew - s) / betas + 0.5 * wtransp / snew
            ysnew = ysprimnew + snew / betas
            srampnew = snew
            ylnew = yl
            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas

        elif wtransp / 2 > vtrimax:
            snew = s + 2 * vtrimax / (ysprim - ylprim)
            ylprimnew = ylprim + (snew - s) / betal
            ysprimnew = ylprimnew
            wrest = wtransp - vtrimax
            wtransp = wrest
            s = snew
            ylprim = ylprimnew
            ysprim = ysprimnew
            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
        # GOTO 1111

    elif ysprim <= ylprim + 0.001:
        snew = ((ys - yl) * s + 2 * wtransp / (1 / betas + 1 / betal)) ** 0.5
        ylprimnew = ylprim
        ysprimnew = ysprim
        ysnew = ysprimnew + snew / betas
        ylnew = ylprimnew - snew / betal
        srampnew = snew
        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
                   
                            
    return snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew
           
#
# snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew=sedbud1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)
#            

# %%

def sedbud0(ys, ysprim, betas, s, yl, ylprim, betal, wtransp):    
# dune evolution fornegative sediment budget
#              print('in sebud0')
#              print('wtransp =', {wtransp})
              if ysprim > ylprim:
                  vtrimax = 0.5 * (ysprim - ylprim) ** 2 * betal
                  if wtransp / 2 <= vtrimax:
                      ylprimnew=ysprim-math.sqrt((ysprim-ylprim)**2-wtransp/betal)
                      snew = s+(ylprimnew-ylprim)*betal
                      ylprim = ylprimnew
                      ylprimnew=ylprim-0.5*wtransp/snew
                      ysprimnew = ysprim
                      ysnew = ys
                      ylnew = yl - 0.5 * wtransp / snew

                  elif wtransp / 2 > vtrimax:
                        snew = s + 2 * vtrimax / (ysprim - ylprim)
                        ylprimnew = ylprim + (snew - s) / betal
                        ysprimnew = ysprim
                        wrest = wtransp - vtrimax
                        wtransp = wrest
                        s = snew
                        snew = math.sqrt(betal * (s * (ysprim - yl) + 2 * wtransp))
                        ysnew = ys
                        ylnew = ylprimnew - snew / betal
                        
              elif ysprim <= ylprim:
                        snew = math.sqrt(betal * (s * (ysprim - yl) + 2 * wtransp))
                        ysnew = ys
                        ylprimnew = ylprim
                        ysprimnew = ysprim
                        ylnew = ylprimnew - snew / betal
                                          
              return snew,ylprimnew,ysprimnew,ylnew,ysnew    
          
# test run
# ys = 147.55
# ysprim=138.8
# ylprim=115.0
# yl = 100.0
# s = 4.14
# betas = 0.47
# betal = 0.26
# wtransp=0.00412004907923848

# snew,ylprimnew,ysprimnew,ylnew,ysnew=sedbud0(ys, ysprim, betas, s, yl, ylprim, betal, wtransp)


# %%

# change of dune volume and geometry due to dune erosion and overwash
def dune(sedbud, vramp, ys, ysprim, betas, swl, r, s, yl, ylprim, betal, dt, qs, ql, qwnet, smax, i, islr, iyear, sramp, vnourd, slr, dfoot, rhov_df, rhov_dt, rhov_db):
    pi = 3.14159


    if r + swl > dfoot:
        erod = (-qs - ql) * dt
#        erod = 0.0  # test to elimate storm erosion
        vrampnew, ysnew, ysprimnew, snew, ylnew, ylprimnew, srampnew = duneeros(vramp, ysprim, betas, betal, s, erod, yl, ylprim, dfoot)
  
                
        if r + swl > dfoot + s:
            deltavback = dt * ql
            ylnew = yl
            ylprimnew = ylprim
            ylnew = yl - deltavback / snew
            ylprim = ylprim - deltavback / snew
          
    elif qwnet > 0:
        deltavramp = dt * qwnet
        vrampnew = vramp + deltavramp

        
        if sedbud == 0:
            vrampmax = (s - 1) ** 2 / (2 * betas)
        else:
            vrampmax = s ** 2 / (2 * betas)

        if vrampnew < vrampmax:
            srampnew = (2 * betas * vrampnew)**0.5
            ysnew = srampnew / betas + ysprim
            ylnew = yl
            ylprimnew = ylprim
            ysprimnew = ysprim
            snew = s
            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
#            wtransp = 0.  # added by Ad
        else:
            wtransp = min(deltavramp, vrampnew - vrampmax)
            vrampnew = max(vrampmax, vramp)
            srampnew = (2 * betas * vrampnew)**0.5
            ysnew = srampnew / betas + ysprim
            ys = ysnew

# morphological evolution depends on sediment budget


            if sedbud == 2:
                if ysprim <= ylprim + 0.001:  # triangular shape
                    snew, ylprimnew, ysprimnew, ylnew, ysnew, vrampnew, srampnew = sedbud2_1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)                 
                elif ysprim > ylprim:    # trapezoidal shape
                    snew, ylprimnew, ysprimnew, ylnew, ysnew, vrampnew, srampnew,wtransp = sedbud2_2(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)                
            elif sedbud == 1:
                snew,ylprimnew,ysprimnew,ylnew,ysnew,vrampnew,srampnew=sedbud1(vramp, ys, ysprim, betas, s, yl, ylprim, betal, smax, wtransp)                 
            elif sedbud == 0:
                snew,ylprimnew,ysprimnew,ylnew,ysnew=sedbud0(ys, ysprim, betas, s, yl, ylprim, betal, wtransp)





                
    else:   # no change in dune profile
        ysnew = ys
        ysprimnew = ysprim
        ylnew = yl
        ylprimnew = ylprim
        snew = s
        srampnew = sramp
        vrampnew = vramp 
        
    if srampnew<0:
        print('BEFORE vegetation')
        print('qwnet', {qwnet})
        print('sedbud', {sedbud})
        print('srampnew', {srampnew})
        print('vrampmax', {vrampmax})
        print('vrampnew', {vrampnew})
        print('wtransp', {wtransp})
        np.pause        

# update vegetation parameters
    rhov_df, rhov_dt, rhov_db = veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, dt,rhov_df,rhov_dt,rhov_db)
        
                
# Update DUNE PARAMETERS for morphological evolution
    ys = ysnew
    ysprim = ysprimnew
    yl = ylnew
    ylprim = ylprimnew
    s = snew
    sramp = srampnew
    vramp = vrampnew 
 

        # ACCOUNT FOR SEA LEVEL RISE
    if i == islr:
        sslr = s - slr * iyear * dt
        yl = yl + slr * iyear * dt / betal
        if ys > ysprim:
            ys = ys - slr * iyear * dt / betas
            if ys < ysprim:
                ys = ysprim
        sramp = (ys - ysprim) * betas
        vramp = 0.5 * (ys - ysprim) * sramp

        s = sslr
        islr = islr + iyear
    
        # ACCOUNT FOR DUNE NOURISHMENT
    if vnourd > 0:
        vrampnew = vramp + vnourd
        vrammax = s**2 / (2 * betas)

        if vrampnew > vrammax:
            wtransp = vrampnew - vrammax
            vrampnew = vrammax
            ysprimnew = ysprim + wtransp / s
            srampnew = s
            ysnew = ysprimnew + s / betas
        else:
            ysprimnew = ysprim
            ysnew = ysprimnew + ((2 * vrampnew) / betas)**0.5
            srampnew = ((2 * vrampnew * betas)**0.5)

        ys = ysnew
        vramp = vrampnew
        ysprim = ysprimnew
        sramp = srampnew
        


        # CALCULATE DUNE VOLUME
    vdune = s * (0.5 * (ylprim - yl) + (ysprim - ylprim)) + (ys - ysprim)**2 * betas / 2
    sramp = srampnew

        
    return vdune,ys,ysprim,vramp,sramp,s,islr,yl,ylprim,rhov_df,rhov_dt,rhov_db
#
#testset up
# sedbud =2
# vramp=13.5
# ys = 145.23
# ysprim = 137.5
# betas = 0.47
# swl = 2.37
# r = 0.0
# s = 3.566
# yl = 100.1
# ylprim = 116.3
# betal = 0.26
# dt = 10800.0
# qs = 9.05685917e-06
# ql = 0.0
# qwnet = 0
# smax = 15
# i = 39786
# islr = 40880
# iyear = 2920
# sramp =3.562
# vnourd=0.0
# slr = 5.8663e-11
# dfoot = 3.0


# dune,ysn,ysprimn,vrampn,srampn,sn,islr,yln,ylprimn,rhov_df,rhov_dt,rhov_db= dune(sedbud, vramp, ys, ysprim, betas, swl, r, s, yl, 
# ylprim, betal, dt, qs, ql, qwnet, smax, i, islr, iyear, sramp, vnourd, slr, dfoot,0,0,0)

# %%
from numba import jit
import numpy as np
# vegetation module to simulate growth/decay
@jit(nopython=True)
def veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, dt,rhov_df,rhov_dt,rhov_db):
    HVmax = 1.0
    HVmin = 0.01 # treshold value
    tveg = 3*24*3600  # duran and moore [2013]
    tveg = 365*24*3600  # needs to be adjusted according to Bonte et al [2021]
    gamv = 10.0  # Duran and Moore[2013] default at 1.0
    # start with ramp
    HV = max((rhov_df)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh = 0.0
    if (abs(ysnew-ys)>0.0):
        dh=(ysnew-ys)*np.sin(betas)
    elif (abs(srampnew-sramp)>0.0):
        dh = 0.5*(srampnew-sramp)
        
#    rhov_dfnew=dt*(1-rhov_df)/tveg-gamv/HV*abs(dh)
    rhov_dfnew=rhov_df+dt*(1-rhov_df)/tveg-gamv*abs(dh)/HVmax
# Bonte
#    rhov_dfnew=rhov_df*1.005-0.005*rhov_df**2
    
    # update
    rhov_df = max(rhov_dfnew,0.0001)
    
    # follow up for the top of the dune
    HV = max((rhov_dt)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh = (snew-s)
        
    rhov_dtnew=rhov_dt+dt*(1-rhov_dt)/tveg-gamv/HVmax*abs(dh)
    # update
    rhov_dt = max(rhov_dtnew,0.0001)
    
    # and the back of the dune
    HV = max((rhov_db)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh=(ylnew-yl)*np.sin(betal)
    rhov_dbnew=rhov_db+dt*(1-rhov_db)/tveg-gamv/HVmax*abs(dh)
    # update
    rhov_db = max(rhov_dbnew,0.0001)
    
    return rhov_df, rhov_dt, rhov_db

# ys = 147.0
# ysnew=147.01
# yl = 100.0
# ylnew = 99.99;
# sramp = 4.0
# srampnew = 4.01
# betas = 0.47
# betal = 0.26
# s = 5
# snew = 5.01
# dt = 10800
# rhov_df = 0.5
# rhov_dt = 0.5
# rhov_db = 0.5

# #
# rhov_df, rhov_dt, rhov_db = veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, dt,rhov_df,rhov_dt,rhov_db)

# %% read CS input file
   
FNAMIN='CSinput.txt'
      
Ins = open(FNAMIN,'r')
line = Ins.readline()
sline=line.split(' ')
NDT =int(sline[0])
DT = float(sline[1])
print(NDT)
print(DT)

line = Ins.readline()
sline=line.split(' ')
TEMP =float(sline[0])
print(TEMP)

line = Ins.readline()
sline=line.split(' ')
SLR =float(sline[0])
print(SLR)

IYEAR=3600*24*365/DT
ISLR=IYEAR    

line = Ins.readline()
sline=line.split(' ')
A =float(sline[0])
print(A)

line = Ins.readline()
sline=line.split(' ')
APRIM =float(sline[0])
print(APRIM)

line = Ins.readline()
sline=line.split(' ')
ASFN =float(sline[0])
print(ASFN)

line = Ins.readline()
sline=line.split(' ')
FNAMWL =sline[0]
print(FNAMWL)

line = Ins.readline()
sline=line.split(' ')
FNAMWIND =sline[0]
print(FNAMWIND)

line = Ins.readline()
sline=line.split(' ')
N =int(sline[0])
print(N)


FNAMP = ["" for j in range(N)]
for j in range(N):
    line = Ins.readline()
    sline=line.split(' ')
    FNAMP[j]=sline[0]
    
Ins.close()

#%% read water level file

Wl = open(FNAMWL,'r')
etal = Wl.readlines()
SWL = np.zeros((NDT));
for j in range(0,NDT):
    SWL[j] = float(etal[j])

Wl.close()
    
#%% read wind conditions file

WIND = open(FNAMWIND,'r')
WSPEEDIN = np.zeros((NDT));
WDIRIN = np.zeros((NDT));
for j in range(0,NDT):
    line = WIND.readline()
    sline=line.split("\t")
    WDIRIN[j] = float(sline[0])
    WSPEEDIN[j] = float(sline[1])
      
WIND.close()



# %% some intermediate plots to check on results

t = np.linspace(0, NDT-1, NDT)*DT/24/3600  # time in days

# Note that even in the OO-style, we use `.pyplot.figure` to create the Figure.
fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(t, SWL, label='tide')  # Plot some data on the axes.
plt.xlabel('t ')
plt.ylabel('WL')
plt.legend()




# %% define variables

IV = np.ones(N)              # COUNTER FOR NOURISHMENTS
QS = np.zeros((NDT,N))       # TRANSPORT FROM DUNE TO BEACH
QL = np.zeros((NDT,N))       # OVERWASH TRANSPORT
MWE = np.zeros((NDT,N))      # POTENTIAL AEOLIAN TRANSPORT RATE
QWNE = np.zeros((NDT,N))     # AEOLIAN TRANSPORT FROM BEACH TO DUNE
QSL = np.zeros((NDT,N))      # TRANSPORT FROM BEACH DUE TO SEA LEVEL RISE
VNOURD = np.zeros((NDT,N))   # NOURISHED VOLUME DUNE
VNOURB = np.zeros((NDT,N))   # NOURISHED VOLUME BEACH
VNOURBAR = np.zeros((NDT,N)) # NOURISHED VOLUME BAR
VDTYEAR =np.zeros((NDT,N))   # CHANGE OF VOLUME IN SEDIMENT BUDGET PREVIOUS YEAR
DB =np.zeros((NDT,N))        # AVERAGE BEACH HEIGHT
R =np.zeros((NDT,N))         # RUNUP
QNOUR = np.zeros((NDT,N))    # ONSHORE TRANSPORT FROM SHOREFACE NOURISHMENTS
YL = np.zeros((NDT,N))       # LANDWARD DUNE FOOT
YSPRIM = np.zeros((NDT,N))   # SEAWARD DUNE CREST END
YS = np.zeros((NDT,N))       # SEAWARD DUNE FOOT
S = np.zeros((NDT,N))        # DUNE HEIGHT
SMAX = np.zeros((N))         # MAX DUNE HEIGHT
BETAL = np.zeros((N))        # LANDWARD SLOPE
BETAS = np.zeros((N))        # TANGENS OF ANGLE OF REPOSE
BETAF = np.zeros((N))        # SWASH SLOPE
D50 = np.zeros((N))          # D50
WFALL = np.zeros((N))        # Fall velocity
SHORN = np.zeros((N))        # SHORE NORMAL
VBEACH = np.zeros((NDT,N))   # Beach volume
VOLBAR = np.zeros((NDT,N))   # bar volume
BUDGET= np.zeros((N))        # initial sediment budget
DCLOS= np.zeros((N))         # depth of closure
DFOOT= np.zeros((N))         # dune foor height
BACTIVE= np.zeros((N))       # active width profile
CIMPACT= np.zeros((N))       # impact coef dune erosion
KW= np.zeros((N))            # aeolian transport coef
CWRATE= np.zeros((N))        # fetch coef
CF= np.zeros((N))            # aeolian fric coef beach
AL= np.zeros((N))            # 
BL= np.zeros((N))            # 
CB= np.zeros((N))            # 
BARNOUR= np.zeros((NDT,N))   # 
VW= np.zeros((NDT,N))        # 
QLSIN= np.zeros((N))         # 
VANR= np.zeros((N))          # 
SEDBUD = np.zeros((NDT,N))   # Sediment budget
YLPRIM = np.zeros((NDT,N))   # 
YG = np.zeros((NDT,N))       # 
VBEACHTOT = np.zeros((NDT,N))# 
VRAMP = np.zeros((NDT,N))    # 
SRAMP = np.zeros((NDT,N))    # 
QWNET = np.zeros((NDT,N))    # 
QSLR = np.zeros((NDT,N))     # 
QLS = np.zeros((NDT,N))      # 
DVOL = np.zeros((NDT,N))     # 
QBAR = np.zeros((NDT,N))     # 
VDUNE = np.zeros((NDT,N))    # 
TETA = np.zeros((NDT,N))     # 
TETAW = np.zeros((NDT,N))    # 
WSPEED = np.zeros((NDT,N))   # 
WDIR = np.zeros((NDT,N))     # 
VDT = np.zeros((NDT,N))      # Volume change in profile
RHOV_DF = np.zeros((NDT,N))  # vegetation density at the front of the dune
RHOV_DT = np.zeros((NDT,N))  # vegetation density on top of the dune
RHOV_DB = np.zeros((NDT,N))  # vegetation density at the back of the dune

#%% read in bed profiles and wave conditions

FNAMWAVE = ["" for j in range(N)] # wave file
HO = np.zeros((NDT,N));
T = np.zeros((NDT,N));
DIR = np.zeros((NDT,N));

VOLAT = np.zeros((N,8))# 
VOLAB = np.zeros((N,8))# 
VOLAD = np.zeros((N,8))# 
VOLRB = np.zeros((N,8))# 
for k in range(N):
    PROF = open(FNAMP[k],'r')
    line = PROF.readline()   # waves file name already known
    sline=line.split(' ')
    FNAMWAVE[k]=sline[0]
    line = PROF.readline()
    sline=line.split()
    YL[0,k] = float(sline[0])
    YSPRIM[0,k] = float(sline[1])
    YS[0,k] = float(sline[2])
    line = PROF.readline()
    sline=line.split()
    S[0,k] = float(sline[0])
    SMAX[k] = float(sline[1])
    line = PROF.readline()
    sline=line.split()
    BETAL[k] = float(sline[0])
    BETAS[k] = float(sline[1])
    BETAF[k] = float(sline[2])
    line = PROF.readline()
    sline=line.split()
    D50[k] = float(sline[0])
    SHORN[k] = float(sline[1])
    line = PROF.readline()
    sline=line.split()
    VBEACH[0,k] = float(sline[0])
    VOLBAR[0,k] = float(sline[1])
    BUDGET[k] = float(sline[2])
    line = PROF.readline()
    sline=line.split()
    DCLOS[k] = float(sline[0])
    DFOOT[k] = float(sline[1])
    BACTIVE[k] = float(sline[2])
    line = PROF.readline()
    sline=line.split()   
    CIMPACT[k] = float(sline[0])
    KW[k] = float(sline[1])
    CWRATE[k] = float(sline[2])
    CF[k] = float(sline[3])
    line = PROF.readline()
    sline=line.split()      
    AL[k] = float(sline[0])
    BL[k] = float(sline[1])
    CB[k] = float(sline[2])
    line = PROF.readline()
    sline=line.split()   
    BARNOUR[0,k] = float(sline[0])
    VW[0,k] = float(sline[1]) 
    line = PROF.readline()
    sline=line.split()   
    QLSIN[k] = float(sline[0])
    line = PROF.readline()
    sline=line.split()   
    VANR[k] = int(sline[0])

                     
    if int(VANR[k])>0:
#         VOLAT = np.zeros((N,int(VANR[k])))# 
         line = PROF.readline()
         sline=line.split()  
         for l in range(int(VANR[k])):    
          VOLAT[k,l] = int(sline[l])
          
#         VOLAD = np.zeros((N,int(VANR[k])))# 
         line = PROF.readline()
         sline=line.split()  
         for l in range(int(VANR[k])):    
          VOLAD[k,l] = float(sline[l])
          
#         VOLAB = np.zeros((N,int(VANR[k])))# 
         line = PROF.readline()
         sline=line.split()  
         for l in range(int(VANR[k])):    
          VOLAB[k,l] = float(sline[l])
          
#         VOLRB = np.zeros((N,int(VANR[k])))# 
         line = PROF.readline()
         sline=line.split()  
         for l in range(int(VANR[k])):    
          VOLRB[k,l] = float(sline[l])

          
    # read in wave conditions      

    WAVE = open(FNAMWAVE[k],'r')
    for j in range(0,NDT):
       line = WAVE.readline()
       sline=line.split()
       HO[j,k] = float(sline[0])
       T[j,k] = float(sline[1])
       DIR[j,k] = float(sline[2])
    WAVE.close()
    
    PROF.close()
    
# %%
# define alongshore transport gradients based on Figure 9 in Hallin et al., 2019
#QLSF =(45.,57.,51.,62.,50.,58.,57.,56.,50.,36.,35.,30.,18.,12.,23.,1.,0.,-10.,-10.,-17.,-18.,-18.,-20.,-20.,-19.,-19.) 
# based on Figure 7
#QLSF =(37.2,44.8,29.8,40.63,33.8,40.5,43.9,42.4,39.,28.4,17.6,24.3,13.8,12.1,17.1,1.,0.,-10.,-10.,-17.,-18.,-18.,-20.,-20.,-19.,-19.) 

# include presence of beach houses
#BH = (1, 1, 1, 1, 0, 0 , 0 , 0 ,0 , 0 , 0 , 0 , 1, 0 , 1, 0 ,0 ,0 ,1,0, 1 , 1, 1 , 1, 1, 1) 


# %% Initialize model variables

I = 0
for K in range(N):
    WFALL[K] = FALVEL(TEMP,D50[K])      
 # Initialize variables
    IV[K] = 1  # COUNTER FOR NOURISHMENTS
    QS[I,K] = 0  # TRANSPORT FROM DUNE TO BEACH
    QL[I,K] = 0  # OVERWASH TRANSPORT
    MWE[I,K]= 0  # POTENTIAL AEOLIAN TRANSPORT RATE
    QWNET[I,K] = 0  # AEOLIAN TRANSPORT FROM BEACH TO DUNE
    QSLR[I,K] = 0  # TRANSPORT FROM BEACH DUE TO SEA LEVEL RISE
    VNOURD[I,K] = 0  # NOURISHED VOLUME DUNE
    VNOURB[I,K] = 0  # NOURISHED VOLUME BEACH
    VNOURBAR[I,K] = 0  # NOURISHED VOLUME BAR
    VDTYEAR[I,K] = 0  # CHANGE OF VOLUME IN SEDIMENT BUDGET PREVIOUS YEAR
    DB[I,K] = 0  # AVERAGE BEACH HEIGHT
    R[I,K] = 0  # RUNUP
    QNOUR[I,K] = 0  # ONSHORE TRANSPORT FROM SHOREFACE NOURISHMENTS
    SEDBUD[I,K]=BUDGET[K]
    YLPRIM[I,K]=YL[I,K]+S[I,K]/BETAL[K]
    if (YLPRIM[I,K]>YSPRIM[I,K]):
        print('PROFILE', K, 'YLPRIM EXCEEDS YSPRIM')
        
    YG[I,K]=YS[I,K]+AL[K]*VBEACH[I,K]+BL[K]   
    DB[I,K]=VBEACH[I,K]/(YG[I,K]-YS[I,K])
    VBEACHTOT[I,K]=YS[I,K]*(DCLOS[K]+DFOOT[K])+(YG[I,K]-YS[I,K])*(DCLOS[K]+DB[I,K])  
    VRAMP[I,K] = 0.5 * (YS[I,K] - YSPRIM[I,K]) ** 2 * BETAS[K] # THIS IS THE MAX
    QLS[I,K] = QLSIN[K]
# Assuming VRAMP is a NumPy array
    if (VRAMP[I,K] > 0):
        SRAMP[I,K] = np.sqrt(2 * BETAS[K] * VRAMP[I,K])
#        if (SRAMP[I,K]>S[0,K]):
#            SRAMP[I,K] = S[0,K]
#            VRAMP[I,K] =  0.5 * (YS[I,K] - YSPRIM[I,K])*(SRAMP[I,K]-DFOOT[K])
    else:
        SRAMP[I,K] = 0
        
#    np.pause
 
# Constants
SQRT = np.sqrt
COS = np.cos
# %% This is the main program
#NDT = 10

# start date of the computations (29-3-96)
Ts = (31*8+30*8+29*8) # number of hours
Tt = 365*8 # total number of time steps in one year (ignoring leap years)
TBH1=4/12  # placement of beach houses
TBH2=10/12  # removal of beach house


for I in range(1,NDT):
    for K in range(N):
        YS[I, K] = YS[I-1, K]
        YL[I, K] = YL[I-1, K]
        S[I, K] = S[I-1, K]
        R[I, K] = R[I-1, K]
        VBEACH[I, K] = VBEACH[I-1, K]
        QS[I, K] = QS[I-1, K]
        YG[I, K] = YG[I-1, K]
        QWNET[I, K] = QWNET[I-1, K]
        QSLR[I, K] = QSLR[I-1, K]
        VNOURB[I, K] = VNOURB[I-1, K]
        VNOURBAR[I, K] = VNOURBAR[I-1, K]
        VW[I, K] = VW[I-1, K]
        QL[I, K] = QL[I-1, K]
        QLS[I, K] = QLS[I-1, K]        
        VOLBAR[I, K] = VOLBAR[I-1, K]
        DVOL[I, K] = DVOL[I-1, K]
        QBAR[I, K] = QBAR[I-1, K]
        DB[I, K] = DB[I-1, K]
        VNOURD[I, K] = VNOURD[I-1, K]
        VRAMP[I, K] = VRAMP[I-1, K]
        YSPRIM[I, K] = YSPRIM[I-1, K]
        YLPRIM[I, K] = YLPRIM[I-1, K]
        SRAMP[I, K] = SRAMP[I-1, K]
        VDUNE[I, K] = VDUNE[I-1, K]
        VBEACHTOT[I, K] = VBEACHTOT[I-1, K]
        BARNOUR[I, K] = BARNOUR[I-1, K]
        QNOUR[I, K] = QNOUR[I-1, K]
        RHOV_DF[I, K] = RHOV_DF[I-1, K]
        RHOV_DT[I, K] = RHOV_DT[I-1, K]
        RHOV_DB[I, K] = RHOV_DB[I-1, K]
 
# CONVERT SIGNIFICANT WAVE HEIGHT TO RMS WAVE HEIGHT
        HO[I, K] = HO[I, K] / 1.414
# MODIFY WAVE HEIGHT WITH RESPECT TO INCIDENT WAVE DIRECTION
        if 90 <= SHORN[K] <= 270:
            UL = SHORN[K] + 90
            LL = SHORN[K] - 90
            if LL <= DIR[I, K] <= UL:
                TETA[I, K] = (DIR[I, K] - SHORN[K]) * 0.01745
            else:
                TETA[I, K] = 0.0
                HO[I, K] = 0.0
                T[I, K] = 0.0
        elif SHORN[K] > 270:
            UL = SHORN[K] + 90 - 360
            LL = SHORN[K] - 90
            if DIR[I, K] >= LL:
                TETA[I, K] = (DIR[I, K] - SHORN[K]) * 0.01745
            elif DIR[I, K] <= UL:
                TETA[I, K] = (DIR[I, K] + 360 - SHORN[K]) * 0.01745
            else:
                TETA[I, K] = 0.0
                HO[I, K] = 0.0
                T[I, K] = 0.0
        elif SHORN[K] < 90:
            UL = SHORN[K] + 90
            LL = SHORN[K] - 90 + 360
            if DIR[I, K] <= UL:
                TETA[I, K] = (DIR[I, K] - SHORN[K]) * 0.01745
            elif DIR[I, K] >= LL:
                TETA[I, K] = (360 - DIR[I, K] + SHORN[K]) * 0.01745
            else:
                TETA[I, K] = 0.0
                HO[I, K] = 0.0
                T[I, K] = 0.0
# Do not allow the waves to become too small
        if T[I, K] <= 0.0:
            T[I, K] = 1.0
        if HO[I, K] <= 0.0:
            HO[I, K] = 0.1
# Calculate wind angle to shore normal, convert to radians, and correct wind speed for non-shoreward angles
        if 80 <= SHORN[K] <= 280:
            UL = SHORN[K] + 80
            LL = SHORN[K] - 80
            if LL <= WDIRIN[I] <= UL:
                TETAW[I, K] = (WDIRIN[I] - SHORN[K]) * 0.01745
                WSPEED[I, K] = WSPEEDIN[I]
            else:
                TETAW[I, K] = 0.0
                WSPEED[I, K] = 0.0
        elif SHORN[K] > 280:
            UL = SHORN[K] + 80 - 360
            LL = SHORN[K] - 80
            if WDIRIN[I] >= LL:
                TETAW[I, K] = (WDIRIN[I] - SHORN[K]) * 0.01745
                WSPEED[I, K] = WSPEEDIN[I]
            elif WDIRIN[I] <= UL:
                TETAW[I, K] = (WDIRIN[I] + 360 - SHORN[K]) * 0.01745
                WSPEED[I, K] = WSPEEDIN[I]
            else:
                TETAW[I, K] = 0.0
                WSPEED[I, K] = 0.0
        elif SHORN[K] < 80:
            UL = SHORN[K] + 80
            LL = SHORN[K] - 80 + 360
            if WDIRIN[I] <= UL:
                TETAW[I, K] = (WDIRIN[I] - SHORN[K]) * 0.01745
                WSPEED[I, K] = WSPEEDIN[I]
            elif WDIRIN[I] >= LL:
                TETAW[I, K] = (360 - WDIRIN[I] + SHORN[K]) * 0.01745
                WSPEED[I, K] = WSPEEDIN[I]
            else:
                TETAW[I, K] = 0.0
                WSPEED[I, K] = 0.0
        WDIR[I, K] = TETAW[I, K]
#        print('I = ', {I})
#
# run some tests
#        HO[I, K] = 0.1  # exclude waves
#        T[I, K] = 1.0
#        QLS[I, K] = (QLSF[K])/365/24/3600. # alongshore transport gradient ( based on figure 9)
#        KW[K] = 2.5 # recalibration
#        QLS[I, K] = 2.0158e-06  # from Caroline
#        SWL[I] = 0.0*SWL[I]   # no tide and surge
#        SLR = 0.0
#        SMAX[K]= 100. # test dynamic max dune height

 
# transport sub aerial part of the profile
        QWNET[I,K],QS[I,K],QL[I,K],VW[I,K],R[I,K] = subar(DT, YS[I,K], YL[I,K], S[I,K], HO[I,K], T[I,K], 
                                           DFOOT[K], D50[K], R[I,K], VBEACH[I,K], CF[K], 
                                           WSPEED[I,K], WDIR[I,K], KW[K], CWRATE[K], CIMPACT[K], 
                                           QS[I,K], SWL[I],BETAF[K], YG[I,K], QWNET[I,K], 
                                           QSLR[I,K], QLS[I,K], VNOURB[I,K], QNOUR[I,K], 
                                           VW[I,K], QL[I,K], K, A, APRIM, ASFN)

# turn of aeolian transport in the presence of a beach house
#        Tfrac = (Ts+I)/(Tt)-math.floor((Ts+I)/(Tt))
#        if (Tfrac>=TBH1 and Tfrac<=TBH2 and BH[K]>0.):  
#            QWNET[I,K]=0.

            
# transport subaqueous part of the profile    
        QBAR[I,K], QNOUR[I,K], QSLR[I,K], DVOL[I,K] = subaq(DT, VOLBAR[I,K], HO[I,K], T[I,K], 
                                           WFALL[K], SLR, BARNOUR[I,K], CB[K], BACTIVE[K])
# specify nourished volumes    
        if int(VANR[K])>0:    
         VNOURD[I,K], VNOURB[I,K], VNOURBAR[I,K], IV[K] = nour(I, IV[K], VOLAT[K,:], VOLAD[K,:], VOLAB[K,:], 
                                                       VOLRB[K,:], VANR[K])
         
# calculate bar evolution    
        VOLBAR[I,K],BARNOUR[I,K] = bar(VOLBAR[I,K], QBAR[I,K], DT, VNOURBAR[I,K], BARNOUR[I,K], QNOUR[I,K])
# create array with volume change in profile


        VDT[I,K]=(-QSLR[I,K]+QLS[I,K])*DT+VNOURD[I,K]+VNOURB[I,K]+VNOURBAR[I,K]
# define sediment budget
        if I < IYEAR:
            SEDBUD[I, K] = BUDGET[K]
        else:
            VDTYEAR[I, K] = sum(VDT[int(I + 1 - IYEAR):I, K]) / IYEAR

        if VDTYEAR[I, K] < 0:
            SEDBUD[I, K] = 0
        elif VDTYEAR[I, K] == 0:
            SEDBUD[I, K] = 1
        else:
            SEDBUD[I, K] = 2
# now compute the dune evolution

        
        VDUNE[I,K],YS[I,K],YSPRIM[I,K],VRAMP[I,K],SRAMP[I,K],S[I,K],ISLR,YL[I,K],YLPRIM[I,K],RHOV_DF[I,K],RHOV_DT[I,K],RHOV_DB[I,K] = dune(
            SEDBUD[I,K], VRAMP[I,K], YS[I,K], YSPRIM[I,K], BETAS[K], 
                          SWL[I], R[I,K], S[I,K], YL[I,K], YLPRIM[I,K], BETAL[K], 
                          DT, QS[I,K], QL[I,K], QWNET[I,K], SMAX[K], I, ISLR, IYEAR, 
                          SRAMP[I,K],VNOURD[I,K], SLR, DFOOT[K], RHOV_DF[I,K],RHOV_DT[I,K],RHOV_DB[I,K])


# compute beach evolution        
        VRAMP[I,K],YS[I,K],YSPRIM[I,K],S[I,K],YL[I,K],YLPRIM[I,K],SRAMP[I,K],VBEACHTOT[I,K],VBEACH[I,K],YG[I,K]=beach(
            QWNET[I,K],QS[I,K],QBAR[I,K],DCLOS[K],DB[I,K],SLR,BACTIVE[K],DT,
            DFOOT[K],VBEACH[I,K],VBEACHTOT[I,K],BETAF[K],YL[I,K],YS[I,K],
            VNOURB[I,K],QLS[I,K],QSLR[I,K],VRAMP[I,K],YSPRIM[I,K],BETAS[K],
            S[I,K],YLPRIM[I,K],BETAL[K],AL[K],BL[K],QNOUR[I,K],SRAMP[I,K])
 
    

# Plot results

# generate figure 14
plt.figure(figsize=(9, 18))
for ip in range(1,K+2):
    plt.subplot(9,3,ip)
    plt.plot(t[1:NDT]/365, VDUNE[1:NDT,ip-1]-VDUNE[1,ip-1])
    plt.plot(t[1:NDT]/365, VBEACHTOT[1:NDT,ip-1]-VBEACHTOT[1,ip-1])
    
# generate figure 14
plt.figure(figsize=(9, 18))
for ip in range(1,K+2):
    plt.subplot(9,3,ip)
    plt.plot(t[1:NDT]/365, VDUNE[1:NDT,ip-1]-VDUNE[1,ip-1])
#    plt.plot(t[1:NDT]/365, VBEACHTOT[1:NDT,ip-1]-VBEACHTOT[1,ip-1])

# plt.subplot(232)
# plt.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K])
# plt.subplot(233)
# plt.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K])
# plt.subplot(234)
# plt.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K])
# plt.subplot(235)
# plt.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K])
# plt.subplot(236)
# plt.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K])


# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# ax.plot(t[1:NDT]/365, VDUNE[1:NDT,K]-VDUNE[1,K], label='change in dune volume')  # Plot some data on the axes.
# ax.plot(t[1:NDT]/365, VBEACHTOT[1:NDT,K]-VBEACHTOT[1,K], label='change in beach volume')  # Plot some data on the axes.
# plt.xlabel('t ')
# plt.ylabel('WS')
# plt.legend()


# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# ax.plot(t, VBEACHTOT[:,K], label='beach volume')  # Plot some data on the axes.
# plt.xlabel('t ')
# plt.ylabel('WS')
# plt.legend()

# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# ax.plot(t, YS[:,K], label='dune foot')  # Plot some data on the axes.
# plt.xlabel('t ')
# plt.ylabel('WS')
# plt.legend()

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(t, QWNET[:,K], label='aeolian transport')  # Plot some data on the axes.
plt.xlabel('t (days)')
plt.ylabel('Q (m^2/s)')
plt.legend()

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(t, RHOV_DF[:,K], label='DUNE RAMP VEGETATION')  # Plot some data on the axes.
plt.xlabel('t (days)')
plt.ylabel('Rveg')
plt.legend()

# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# ax.plot(t, WDIR[:,K], label='WIND SPEED')  # Plot some data on the axes.
# plt.xlabel('t ')
# plt.ylabel('WSPEED')
# plt.legend()

# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
Y = (0, YL[I,K], YLPRIM[I,K], YSPRIM[I,K], YSPRIM[I,K]+1.0, YS[I,K], YG[I,K],YG[I,K]+500.0)
Z = (DFOOT[K], DFOOT[K] , S[I,K]+DFOOT[K], S[I,K]+DFOOT[K],SRAMP[I,K]+DFOOT[K],DFOOT[K], 0.0 ,-DCLOS[K])
# ax.plot(Y, Z, label='PROFILE')  # Plot some data on the axes.
# plt.xlabel('Y (M) ')
# plt.ylabel('Z (M)')
# plt.legend()

# plot difference
# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# Y0 = (0, YL[0,K], YLPRIM[0,K], YSPRIM[0,K], YSPRIM[0,K]+1.0, YS[0,K], YG[0,K],YG[0,K]+500.0)
# Z0 = (DFOOT[K], DFOOT[K] , S[0,K], S[0,K],SRAMP[0,K],DFOOT[K], 0.0 ,-DCLOS[K])
# ax.plot(Y, Z, label='Final profile')  # Plot some data on the axes.
# ax.plot(Y0,Z0, label='Initial Profile') 
# plt.xlabel('Y (M) ')
# plt.ylabel('Z (M)')
# plt.legend()


fig, ax = plt.subplots(figsize=(6, 4), layout='constrained')
Y0 = (0, YL[0,K], YLPRIM[0,K], YSPRIM[0,K], YSPRIM[0,K]+1.0, YS[0,K], YG[0,K],YG[0,K]+500.0)
Z0 = (DFOOT[K], DFOOT[K] , S[0,K]+DFOOT[K], S[0,K]+DFOOT[K],SRAMP[0,K]+DFOOT[K],DFOOT[K], 0.0 ,-DCLOS[K])
# add vegetation
YV = (0, YL[I,K], YLPRIM[I,K], YSPRIM[I,K], YSPRIM[I,K]+1.0, YS[I,K])
VF = (DFOOT[K]+np.sqrt(RHOV_DB[I,K]), DFOOT[K] +np.sqrt(RHOV_DB[I,K]), S[I,K] +DFOOT[K]+np.sqrt(RHOV_DT[I,K]), S[I,K] +DFOOT[K]
      +np.sqrt(RHOV_DT[I,K]),SRAMP[I,K] +DFOOT[K]+np.sqrt(RHOV_DF[I,K]),DFOOT[K])
ax.plot(Y, Z, label='Final profile')  # Plot some data on the axes.
ax.plot(Y0,Z0, label='Initial Profile') 
ax.plot(YV,VF, label='vegetation') 
plt.xlabel('Y (M) ')
plt.ylabel('Z (M)')
plt.legend()
axes = plt.gca()
axes.set_xlim([50,500])
axes.set_ylim([-2,20])


# fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# ax.plot(t, VOLBAR[:,K], label='bar volume')  # Plot some data on the axes.
# plt.xlabel('t ')
# plt.ylabel('WS')
# plt.legend()


from matplotlib import cm
from matplotlib.colors import LightSource
# create a spatial plot

# interpolate y-axis to regular grid
YI = np.linspace(0, 4000, 500)
ny = np.size(YI)
XI = np.linspace(0,25, 25)
nx = np.size(XI)
ZI0 = np.zeros((nx,ny))
ZI1 = np.zeros((nx,ny))
# add vegetation
YV_DF = np.zeros((nx))
ZV_DF =  np.zeros((nx))
TVI =  np.zeros((nx))
ip = 0

# create a spatial plot
# add shore line position
# YS0=(2727,2272,2090,1636,1363,1272,1090,909,818,727,590,545,363,363,363,363,363,363,363,363,363,363,363,363,363,363)


# for i in range(25):

#     Y0 = (0, YS0[i], YL[0,i]+YS0[i], YLPRIM[0,i]+YS0[i], YSPRIM[0,i]+YS0[i], 
#           YSPRIM[0,i]+1.0+YS0[i], YS[0,i]+YS0[i], YG[0,i]+YS0[i],YG[0,i]+500.0+YS0[i])
#     Z0 = (DFOOT[i],DFOOT[i], DFOOT[i] , S[0,i], S[0,i],SRAMP[0,i],DFOOT[i], 0.0 ,-DCLOS[i])
#     ZI0[ip,:] = np.interp(YI,Y0,Z0)
#     Y1 = (0,YS0[i], YL[I,i]+YS0[i], YLPRIM[I,i]+YS0[i], YSPRIM[I,i]+YS0[i],
#           YSPRIM[I,i]+1.0+YS0[i], YS[I,i]+YS0[i], YG[I,i]+YS0[i],YG[I,i]+500.0+YS0[i])
#     Z1 = (DFOOT[i],DFOOT[i], DFOOT[i] , S[I,i], S[I,i],SRAMP[I,i],DFOOT[i], 0.0 ,-DCLOS[i])
#     ZI1[ip,:] = np.interp(YI,Y1,Z1)

#     YV_DF[ip]=YSPRIM[I,i]
#     ZV_DF[ip]=np.sqrt(RHOV_DF[I,i])+S[I,i]
#     TVI[ip] = ip
#     ip = ip +1

# YI, XI = np.meshgrid(YI, XI) 
   
# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))    
# surf = ax.plot_surface(YI, XI, ZI0, cmap=cm.viridis,
#                        linewidth=0, antialiased=False) 
# surf.set_clim(vmin=-3, vmax=10)
# ax.set_zlim((-5, 20))

# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))    
# surf = ax.plot_surface(YI, XI, ZI1, cmap=cm.viridis,
#                        linewidth=0, antialiased=False) 
# surf.set_clim(vmin=-3, vmax=10)
# ax.set_zlim((-5, 20))


# plt.show()

# # create a temporal plot


# # interpolate y-axis to regular grid
YI = np.linspace(0, 500, 500)
ny = np.size(YI)
TI = np.linspace(0,1000, 1000)
nt = np.size(TI)
ZI = np.zeros((nt,ny))
# add vegetation
YV_DF = np.zeros((nt))
ZV_DF =  np.zeros((nt))
TVI =  np.zeros((nt))
ip = 0


for i in range(0,60000,60):

    Y0 = (0, YL[i,K], YLPRIM[i,K], YSPRIM[i,K], YSPRIM[i,K]+1.0, YS[i,K], YG[i,K],YG[i,K]+500.0)
    Z0 = (DFOOT[K], DFOOT[K] , S[i,K], S[i,K],SRAMP[i,K],DFOOT[K], 0.0 ,-DCLOS[K])
    ZI[ip,:] = np.interp(YI,Y0,Z0)
    YV_DF[ip]=YSPRIM[i,K]
    ZV_DF[ip]=np.sqrt(RHOV_DF[i,K])+S[i,K]
    TVI[ip] = ip
    ip = ip +1

YI, TI = np.meshgrid(YI, TI)    
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))    
surf = ax.plot_surface(YI, TI, ZI, cmap=cm.viridis,
                        linewidth=0, antialiased=False) 
surf.set_clim(vmin=-3, vmax=10)

ax.scatter(YV_DF,TVI, zs=7,zdir='z', c=ZV_DF, cmap='YlGn')
ax.stem(YV_DF[0:1000:10],TVI[0:1000:10],ZV_DF[0:1000:10],linefmt='none',bottom=DFOOT[K])
ax.stem([200,200],[100, 800],[10, 8],linefmt='green')
# Customize the z axis.
ax.set_zlim((-5, 20))

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(YI[1,:], ZI[1,:],YI[100,:], ZI[100,:],YI[200,:], ZI[200,:],YI[500,:], ZI[500,:],YI[999,:], ZI[999,:] )  # Plot some data on the axes.
#ax.plot(YI[1:100:1000,:], ZI[1:100:1000,:] )  # Plot some data on the axes.
plt.xlabel('y (m) ')
plt.ylabel('z (m)')
plt.legend()

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(t, S[:,K], label='Dune Height')  # Plot some data on the axes.
plt.xlabel('t (days)')
plt.ylabel('S (m)')
plt.legend()

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(t, YG[:,K]-YS[:,K], label='Lveg')  # Plot some data on the axes.
plt.xlabel('t ')
plt.ylabel('Lveg (m)')
plt.legend()


# add vegetation



# pc = axs.scatter(Y0,X0 , c=Z0, cmap='RdBu_r')
# pc = axs.scatter(Y0,X0 , c=Z0, cmap='cividis')
# XV = np.ones(np.size(YV))*i
# YV = (0, YL[i,K], YLPRIM[i,K], YSPRIM[i,K], YSPRIM[i,K]+1.0, YS[i,K])
# ZV = (DFOOT[K]+np.sqrt(RHOV_DB[i,K]), DFOOT[K] +np.sqrt(RHOV_DB[i,K]), S[0,K] +np.sqrt(RHOV_DT[i,K]), S[0,K] 
#           +np.sqrt(RHOV_DT[i,K]),SRAMP[0,K] +np.sqrt(RHOV_DF[i,K]),DFOOT[K])
# pc = axs.scatter(YV,XV , c=ZV, cmap='YlGn')
# plt.clim(-15, 10)    
# axes = plt.gca()
# axes.set_xlim([50,300]) 
# axes.set_ylim([0,60000])
  
# fig.colorbar(pc, ax=axs[0, 0])
# axs[0, 0].set_title('pcolormesh()')


# create an animation

