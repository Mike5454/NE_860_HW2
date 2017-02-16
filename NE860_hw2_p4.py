from __future__ import division
import matplotlib.pyplot as plt
import re
import numpy as np

#########################################################
#           Constants & Cross sections
#########################################################

def sigma(name):
    energy=[]
    sigma=[]
    Eng=np.logspace(0,2,2000)
#    Eng=np.linspace(1,100,2000)
    F = open(name, 'r')
    line = F.readline().strip()
    for line in F:
        column = line.split()
        energy.append(column[0])
        sigma.append(float(column[1])*1E-24)
    data = np.interp(Eng, energy, sigma)
    return data
    
U238t= sigma('U238_sigma_t.txt')
U238es= sigma('U238_sigma_es.txt')
H1t= sigma('H1_sigma_t.txt')
H1es= sigma('H1_sigma_es.txt')

###############################################################################

def scatter_probability(E_i, E_j, alpha) :
    p = (1.0/E_j/(1.0-alpha))*1.0*((E_i >= alpha*E_j))
    return p
    
def compute_spectrum(E, Sigma_t, Sigma_s, alpha) :
    N_E = len(E)
    phi = np.zeros(N_E)
    phi[N_E-1] = 1.0/(Sigma_t[N_E-1]-Sigma_s[N_E-1])
    for i in range(N_E-2, -1, -1) :
        Q_i = 0.0
        for j in range(N_E-1, i, -1) :
            dE = E[j] - E[j-1]
            E_bar = np.sqrt(E[j]*E[j-1])
            Q_i += phi[j]*dE*(NH1*H1es[j] * scatter_probability(E[i], E_bar, alpha[0]) + NU238*U238es[j]*scatter_probability(E[i], E_bar, alpha[1]))
        phi[i] = Q_i / (Sigma_t[i])
#        print i, phi[i]
    return phi

#########################################################
#                   Inputs
#########################################################

ratio = [1,10,100,1000,10000,100000]
# This will give us the 1000 to 1 ratio for hydrogen to Uranium for an arbitrary density
for q in ratio:
    NH1 = 8.988E-5* 6.022E23 # (rho * Na)/M including conversion
            
    NU238 = NH1/q # (rho * Na)/M including conversion
    
    A = [1,238]
    alpha =[]
    for i in A:
        alpha.append(((i-1)/(i+1))**2)
    sigmat=[]
    sigmael=[]
    sigma_srE=[]
    sigma_sr = 0
    sigma_sk= 0
    E = np.logspace(0, 2, 2000)
    for i in range(0,len(E)):
        sigmat.append(NH1*H1t[i] + NU238*U238t[i])
        sigmael.append(NH1*H1es[i] + NU238*U238es[i])
        sigma_sr = U238es[i]*NU238 + sigma_sr
        sigma_sk = H1es[i]*NH1 + sigma_sk
        sigma_srE.append(U238es[i]*NU238)
    #print sigmael
    #print sigmat 
    phi = compute_spectrum(E, sigmat, sigmael, alpha)
    ############################################################
    #           Narrow and Wide calculations
    ############################################################
    phin=[]
    phiwr=[]
    for i in range(0,len(E)):
        phin.append(((sigma_sr/100)+sigma_sk)/(sigmat[i]*E[i]))
        phiwr.append(sigma_sk/((sigmat[i]-sigma_srE[i])*E[i]))
    
    ###############################################################################
    #              Calculate sigma U238 and Backgrount Cross-section
    ###############################################################################
    # Find the group capture cross section of U238 for the numerical solution
    i=0
    num=0
    nnum=0
    denom=0
    num2=0
    nnum2=0
    denom2=0
    num3=0
    nnum3=0
    denom3=0
    if q==1:
        sigmab=[]
        sigmaU238=[]
    while E[i] <= 12:
        num = ((U238t[i]-U238es[i])*1E24 *phi[i]*(E[i+1]-E[i])) + num
        denom = (phi[i]*(E[i+1]-E[i])) + denom
        nnum = NH1*H1t[i]*1E24/(E[i+1]-E[i]) + nnum
        i= i +1
    sigmaU238.append(num/denom)
    sigmab.append(nnum/(NU238))
    while E[i]>=12 and E[i] <= 28:
        num2 = ((U238t[i]-U238es[i])*1E24 *phi[i]*(E[i+1]-E[i])) + num2
        denom2 = (phi[i]*(E[i+1]-E[i])) + denom2
        nnum2 = NH1*H1t[i]*1E24/(E[i+1]-E[i]) + nnum2
        i= i +1
    sigmaU238.append(num2/denom2)
    sigmab.append(nnum2/(NU238))
    while E[i]>=28 and E[i] <= 50:
        num3 = ((U238t[i]-U238es[i])*1E24 *phi[i]*(E[i+1]-E[i])) + num3
        denom3 = (phi[i]*(E[i+1]-E[i])) + denom3
        nnum3 = NH1*H1t[i]*1E24/(E[i+1]-E[i]) + nnum3
        i= i +1
    sigmaU238.append(num3/denom3)
    sigmab.append(nnum3/(NU238))
    # Repeat process for Narrow Resonance
    i=0
    num=0
    denom=0
    num2=0
    denom2=0
    num3=0
    denom3=0
    if q==1:
        sigmaU238n=[]    
    while E[i] <= 12:
        num = ((U238t[i]-U238es[i])*1E24 *phin[i]*(E[i+1]-E[i])) + num
        denom = (phin[i]*(E[i+1]-E[i])) + denom
        i= i +1
    sigmaU238n.append(num/denom)
    while E[i]>=12 and E[i] <= 28:
        num2 = ((U238t[i]-U238es[i])*1E24 *phin[i]*(E[i+1]-E[i])) + num2
        denom2 = (phin[i]*(E[i+1]-E[i])) + denom2
        i= i +1
    sigmaU238n.append(num2/denom2)
    while E[i]>=28 and E[i] <= 50:
        num3 = ((U238t[i]-U238es[i])*1E24 *phin[i]*(E[i+1]-E[i])) + num3
        denom3 = (phin[i]*(E[i+1]-E[i])) + denom3
        i= i +1
    sigmaU238n.append(num3/denom3)
    # Repeat process for the Wide Resonance
    i=0
    num=0
    denom=0
    num2=0
    denom2=0
    num3=0
    denom3=0
    if q==1:
        sigmaU238wr=[]    
    while E[i] <= 12:
        num = ((U238t[i]-U238es[i])*1E24 *phiwr[i]*(E[i+1]-E[i])) + num
        denom = (phiwr[i]*(E[i+1]-E[i])) + denom
        i= i +1
    sigmaU238wr.append(num/denom)
    while E[i]>=12 and E[i] <= 28:
        num2 = ((U238t[i]-U238es[i])*1E24 *phiwr[i]*(E[i+1]-E[i])) + num2
        denom2 = (phiwr[i]*(E[i+1]-E[i])) + denom2
        i= i +1
    sigmaU238wr.append(num2/denom2)
    while E[i]>=28 and E[i] <= 50:
        num3 = ((U238t[i]-U238es[i])*1E24 *phiwr[i]*(E[i+1]-E[i])) + num3
        denom3 = (phiwr[i]*(E[i+1]-E[i])) + denom3
        i= i +1
    sigmaU238wr.append(num3/denom3)
    
    
#G = open("sigmacU238.txt", 'w')
#for i in sigmaU238:
#    G.write(str(i)+'\n')
#H = open("sigmacU238n.txt", 'w')
#for i in sigmaU238n:
#    H.write(str(i)+'\n')
#I = open("sigmacU238wr.txt", 'w')
#for i in sigmaU238wr:
#    I.write(str(i)+'\n')
#j = open("sigmab.txt", 'w')
#for i in sigmab:
#    j.write(str(i)+'\n')
    
print " Relative errors are given for the 3 energy groups of the first ratio"
print "followed by the errors of the 3 energy groups for the second ratio, and so on"
print "This is done for the narrow results and then the wide results."
rel_nar = []
for i in range(0,18):
    er = (sigmaU238n[i]-sigmaU238[i])/sigmaU238[i]
    rel_nar.append(er)
    
rel_wide = []
for i in range(0,18):
    err = (sigmaU238wr[i]-sigmaU238[i])/sigmaU238[i]
    rel_wide.append(err)
print "\n"
print "Relative errors for Narrow"
print rel_nar[0:3]
print rel_nar[3:6]
print rel_nar[6:9]
print rel_nar[9:12]
print rel_nar[12:15]
print rel_nar[15:18]
print "Relative errors for Wide"
print rel_wide[0:3]
print rel_wide[3:6]
print rel_wide[6:9]
print rel_wide[9:12]
print rel_wide[12:15]
print rel_wide[15:18]

