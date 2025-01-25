# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:01:39 2025

@author: assensme
"""

#Plots used for the article, per scenario
    
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import numpy as np
import matplotlib.pyplot as plt
import PEC_cell
voltage_range=np.arange(-0.0,-1.7,-0.05)
j_range=np.zeros(len(voltage_range))
FE_range=np.arange(0,100)

for i in range(len(voltage_range)):
    P=PEC_cell.PEC_Cell()
    P.Bandgap_message=False
    P.Eg_completed=False
    P.message_spectrum_loaded=False
    P.plot_curves=False
    P.set_voltage_max_ethylene=voltage_range[i]
    P.set_max_ethylene=61.2
    P.Calculate_setup(scenario='VariableCO2RRcatalyst')
    j_range[i]=np.max(P.Sparse_Coverage_J_matrix)
    print(j_range[i])

plt.figure()
plt.plot(voltage_range,j_range)
plt.xlabel(r'$\eta_{\rm CO2RR}$ (V vs RHE)')
plt.ylabel(r'Maximum j$_{\rm s}$ (mA/cm$^2$)')




X=np.meshgrid(voltage_range,FE_range)[0]
Xj=np.meshgrid(j_range,FE_range)[0]
Y=np.meshgrid(voltage_range,FE_range)[1]
umol_production=np.zeros(np.shape(X))
STE_production=np.zeros(np.shape(X))


for i in range(np.shape(X)[1]):
    for j in range(np.shape(X)[0]):
        umol_production[j,i]=P.mol_ethylene(j_range[i],FE_range[j])
        STE_production[j,i]=P.STE(j_range[i],FE_range[j])
        
plt.figure()

unit=r' $\mu$mol/h/cm$^2$ C$_2$H$_4$'
cmap='Purples_r'
levels=np.arange(0,80,10)
clabel=r'C$_2$H$_4$ production rate ($\mu$mol/h/cm$^2$)'
figcontour, ax = plt.subplots()
  
    
CS1 = ax.contourf(X, Y, umol_production,levels=np.arange(0,80,1),cmap=cmap)
CS2 = ax.contour(X, Y, umol_production,levels=levels,colors='w')

cbar=figcontour.colorbar(CS1,ticks=np.arange(0,80,10))    
cbar.ax.set_ylabel(clabel)

plt.xlabel(r'$\eta_{\rm CO2RR}$ (V vs RHE)')

plt.ylabel('Faradaic Efficiency (%)')

plt.plot(-0.61,61.2,'ok')
plt.plot(-0.61,70,'*k')
plt.plot(-0.4,61.2,'*k')

plt.plot(-1.55,60.7,'Xk')
plt.plot(-0.68,57,'^k')
plt.plot(-0.70,59.5,'Pk')
plt.plot(-1.05,51.5,'sk')
plt.plot(-1.05,38.7,'dk')
plt.plot(-1.00,33.2,'pk')
plt.figure()

boxtext='STE Efficiency: '
unit='%'
levels=np.arange(0,30,5)
cmap='Greens_r'
clabel='STE efficiency (%)'
figcontour, ax = plt.subplots()
    
CS1 = ax.contourf(X, Y, STE_production,levels=np.arange(0,30,0.5),cmap=cmap)
CS2 = ax.contour(X, Y, STE_production,levels=levels,colors='w')

cbar=figcontour.colorbar(CS1,ticks=levels)    
cbar.ax.set_ylabel(clabel)

plt.xlabel(r'$\eta_{\rm CO2RR}$ (V vs RHE)')

plt.ylabel('Faradaic Efficiency (%)')

plt.plot(-0.61,61.2,'ok')
plt.plot(-0.61,70,'*k')
plt.plot(-0.4,61.2,'*k')

plt.plot(-1.55,60.7,'Xk')
plt.plot(-0.68,57,'^k')
plt.plot(-0.70,59.5,'Pk')
plt.plot(-1.05,51.5,'sk')
plt.plot(-1.05,38.7,'dk')
plt.plot(-1.00,33.2,'pk')
plt.show()

print('Isocurves')