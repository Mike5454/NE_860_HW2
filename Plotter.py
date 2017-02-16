import numpy as np
import matplotlib.pyplot as plt

A = open('results10.txt','r')
B = open('results100.txt','r')
C = open('results1000.txt','r')
D= open('NR10.txt', 'r')
E= open('NR100.txt', 'r')
F= open('NR1000.txt', 'r')
G= open('Wr10.txt', 'r')
H= open('Wr100.txt', 'r')
I= open('Wr1000.txt', 'r')
J=open('sigmacU238.txt', 'r')
K=open('sigmacU238n.txt', 'r')
L=open('sigmacU238wr.txt', 'r')
M=open('sigmab.txt', 'r')

Results10=[]
for i in A:
    Results10.append(i)
Results100=[]
for i in B:
    Results100.append(i)
Results1000=[]
for i in C:
    Results1000.append(i)
    
NR10=[]
for i in D:
    NR10.append(i)
NR100=[]
for i in E:
    NR100.append(i)
NR1000=[]
for i in F:
    NR1000.append(i)
    
WR10=[]
for i in G:
    WR10.append(i)
WR100=[]
for i in H:
    WR100.append(i)
WR1000=[]
for i in I:
    WR1000.append(i)
    
Eng = np.logspace(0,2,5000)

sigmaU238=[]
for i in J:
    sigmaU238.append(i)
    
sigmaU238n=[]
for i in K:
    sigmaU238n.append(i)
    
sigmaU238wr=[]
for i in L:
    sigmaU238wr.append(i)
sigmab=[]
for i in M:
    sigmab.append(i)


#########################################################################################
#                   Comparisons Ratio
#########################################################################################
fig = plt.Figure
plt.figure(figsize = (15,10.5))

plt.xlabel("Energy (eV)",  fontname="Arial", fontsize=30)
plt.xscale('log')
#plt.xlim(9, 10)

plt.tick_params(which='major', length=10, labelsize=25)
plt.tick_params(which='minor', length=7)
#plt.ticklabel_format(ScalarFormatter(useOffset=False))

plt.ylabel(r'$\phi$ (cm$^{-2}$s$^{-1}$ev$^{-1})$', fontname="Arial", fontsize=30)
plt.yscale('log')
#plt.ylim(6, 11)

plt.grid(True, which='minor', color='lightgrey', linestyle='-')
plt.grid(True, which='major', color='dimgrey', linestyle='-')

plt.title ("Neutron Spectrum",fontsize=30)
plt.rc('font',family='Arial')

p1 = plt.plot(Eng, Results10, 'k-', label = r'Ratio $\frac{10H}{1U}$', linewidth = 4)
p2 = plt.plot(Eng, Results100, 'b-', label = r'Ratio $\frac{100H}{1U}$', linewidth = 4)
p3 = plt.plot(Eng, Results1000, 'r-', label = r'Ratio $\frac{1000H}{1U}$', linewidth = 4)

plt.legend(loc=1,prop={'size':20})
plt.show()



#########################################################################################
#                       Comparisons NR & WR
#########################################################################################
fig = plt.Figure
plt.figure(figsize = (15,10.5))

plt.xlabel("Energy (eV)",  fontname="Arial", fontsize=30)
plt.xscale('log')
#plt.xlim(9, 10)

plt.tick_params(which='major', length=10, labelsize=25)
plt.tick_params(which='minor', length=7)
#plt.ticklabel_format(ScalarFormatter(useOffset=False))

plt.ylabel(r'$\phi$ (cm$^{-2}$s$^{-1}$ev$^{-1})$', fontname="Arial", fontsize=30)
plt.yscale('log')
#plt.ylim(6, 11)

plt.grid(True, which='minor', color='lightgrey', linestyle='-')
plt.grid(True, which='major', color='dimgrey', linestyle='-')

plt.title ("Neutron Spectrum",fontsize=30)
plt.rc('font',family='Arial')

p1 = plt.plot(Eng, Results10, 'k-', label = r'Ratio $\frac{10H}{1U}$', linewidth = 4)
p2 = plt.plot(Eng, NR10, 'b-', label = "Narrow Resonance", linewidth = 4)
p3 = plt.plot(Eng, WR10, 'r-', label = "Wide Resonance", linewidth = 4)

plt.legend(loc=1,prop={'size':20})
plt.show()


#########################################################################################
#               Extract sigma data for all ratios
#########################################################################################

Uc1=[]
Uc10=[]
Uc100=[]
Uc1000=[]
Uc10000=[]
Uc100000=[]
Ub1=[]
Ub10=[]
Ub100=[]
Ub1000=[]
Ub10000=[]
Ub100000=[]
for i in range(0,3):
    Uc1.append(sigmaU238[i])
    Ub1.append(sigmab[i])
for i in range(3,6):
    Uc10.append(sigmaU238[i])
    Ub10.append(sigmab[i])
for i in range(6,9):
    Uc100.append(sigmaU238[i])
    Ub100.append(sigmab[i])
for i in range(9,12):
    Uc1000.append(sigmaU238[i])
    Ub1000.append(sigmab[i])
for i in range(12,15):
    Uc10000.append(sigmaU238[i])
    Ub10000.append(sigmab[i])
for i in range(15,18):
    Uc100000.append(sigmaU238[i])
    Ub100000.append(sigmab[i])

#########################################################################################
#       Plotting Sigma capture of U235 vs sigma background
#########################################################################################
fig = plt.Figure
plt.figure(figsize = (15,10.5))

plt.xlabel(r'$\sigma_b$',  fontname="Arial", fontsize=30)
plt.xscale('log')
#plt.xlim(9, 10)

plt.tick_params(which='major', length=10, labelsize=25)
plt.tick_params(which='minor', length=7)
#plt.ticklabel_format(ScalarFormatter(useOffset=False))

plt.ylabel(r'$\sigma^{U238}_c$', fontname="Arial", fontsize=30)
plt.yscale('log')
#plt.ylim(6, 11)

plt.grid(True, which='minor', color='lightgrey', linestyle='-')
plt.grid(True, which='major', color='dimgrey', linestyle='-')

plt.title ("Capture Cross Section",fontsize=30)
plt.rc('font',family='Arial')

p1 = plt.step(Ub1, Uc1, 'k-', label = r'$\frac{1H}{1U}$', linewidth = 4)
p2 = plt.step(Ub10, Uc10, 'm-', label = r'$\frac{10H}{1U}$', linewidth = 4)
p3 = plt.step(Ub100, Uc100, 'r-', label = r'$\frac{100H}{1U}$', linewidth = 4)
p4 = plt.step(Ub1000, Uc1000, 'c-', label = r'$\frac{1000H}{1U}$', linewidth = 4)
p5 = plt.step(Ub10000, Uc10000, 'g-', label = r'$\frac{10000H}{1U}$', linewidth = 4)
p6 = plt.step(Ub100000, Uc100000, 'b-', label = r'$\frac{100000H}{1U}$', linewidth = 4)

plt.legend(loc=4,prop={'size':20})
plt.show()
