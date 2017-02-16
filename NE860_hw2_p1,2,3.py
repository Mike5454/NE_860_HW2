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

ratio = 10
# This will give us the 1000 to 1 ratio for hydrogen to Uranium for an arbitrary density
NH1 = 8.988E-5* 6.022E23 # (rho * Na)/M including conversion
        
NU238 = NH1/ratio # (rho * Na)/M including conversion

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
    
############# Normalizing ##############
#phi_norm=[]
#phin_norm=[]
#phiwr_norm=[]
#for i in phi:
#    phi_norm.append(i/phi[0])
#for i in phin:
#    phin_norm.append(i/phin[0])
#for i in phiwr:
#    phiwr_norm.append(i/phiwr[0])

###############################################################################
#              Calculate sigma U238 and Backgrount Cross-section
###############################################################################
i=0
val=0
sigmaU238=[]
while E[i] <= 12:
    val = ((U238t[i]-U238es[i])*10**24 *phi[i]*(E[i]-E[i+1])/(phi[i]*(E[i]-E[i+1]))) + val
    i= i +1
    print i
sigmaU238.append(val)
while E[i]>=12 and E[i] <= 28:
    val = ((U238t[i]-U238es[i])*10**24 *phi[i]*(E[i]-E[i+1])/(phi[i]*(E[i]-E[i+1]))) + val
    i= i+1
sigmaU238.append(val)
while E[i]>=28 and E[i] <= 50:
    val = ((U238t[i]-U238es[i])*10**24 *phi[i]*(E[i]-E[i+1])/(phi[i]*(E[i]-E[i+1]))) + val
    i= 1+1
sigmaU238.append(val)
##################################################
# Writing to text file to be graphed in plotter    
##################################################
#G = open("Wr"+str(ratio)+".txt", 'w')
#for i in phiwr:
#    G.write(str(i)+'\n')
#H = open("NR"+str(ratio)+".txt", 'w')
#for i in phin:
#    H.write(str(i)+'\n')
#I = open("Results"+str(ratio)+".txt", 'w')
#for i in phi:
#    I.write(str(i)+'\n')
###############################################################################
#                            Test Plots 
###############################################################################
#fig = plt.Figure
#plt.figure(figsize = (15,10.5))
#
#plt.xlabel("Energy (eV)",  fontname="Arial", fontsize=30)
#plt.xscale('log')
##plt.xlim(9, 10)
#
#plt.tick_params(which='major', length=10, labelsize=25)
#plt.tick_params(which='minor', length=7)
##plt.ticklabel_format(ScalarFormatter(useOffset=False))
#
#plt.ylabel(r'$\phi$ (cm$^{-2}$s$^{-1}$ev$^{-1})$', fontname="Arial", fontsize=30)
#plt.yscale('log')
##plt.ylim(6, 11)
#
#plt.grid(True, which='minor', color='lightgrey', linestyle='-')
#plt.grid(True, which='major', color='dimgrey', linestyle='-')
#
#plt.title ("Neutron Spectrum",fontsize=30)
#plt.rc('font',family='Arial')
#
#p1 = plt.plot(E, phi, 'k-', label = r'Ratio $\frac{10H}{1U}$', linewidth = 4)
#
#plt.legend(loc=2,prop={'size':20})
#plt.show()



#
#fig = plt.Figure
#plt.figure(figsize = (15,10.5))
#plt.text(7950, -0.05, r'8x10$^3$',  fontsize = 25)
#
#plt.xlabel("Energy (eV)",  fontname="Arial", fontsize=30)
#plt.xscale('log')
#plt.xlim(8E3, 1E4)
#
#plt.tick_params(which='major', length=10, labelsize=25)
#plt.tick_params(which='minor', length=7)
##plt.ticklabel_format(ScalarFormatter(useOffset=False))
#
#plt.ylabel(r'$\phi$ (cm$^{-2}$s$^{-1}$ev$^{-1})$', fontname="Arial", fontsize=30)
## plt.yscale('log')
##plt.ylim(6, 11)
#
#plt.grid(True, which='minor', color='lightgrey', linestyle='-')
#plt.grid(True, which='major', color='dimgrey', linestyle='-')
#
#plt.title ("Neutron Spectrum",fontsize=30)
#plt.rc('font',family='Arial')
#
#p1 = plt.plot(E, phi, 'k-', label = 'U-238 Spectrum', linewidth = 4)
#
#plt.legend(loc=2,prop={'size':20})
#plt.show()
