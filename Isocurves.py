# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:01:39 2025

This script calculates js using different overpotentials of the CO2RR catalyst under AM1.5G using tandem solar cells. 
Together with the faradaic efficiency, which is varied from 0 to 100, isocurves of production rate can be make, generating Figure 5 and Figure S7. 

The OER catalyst is assumed to be NiFeOx, R_solution=5 

@author: assensme
"""

    
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import numpy as np
import matplotlib.pyplot as plt
import PEC_cell
voltage_range=np.arange(-0.0,-1.7,-0.05)
j_range=np.zeros(len(voltage_range))
FE_range=np.arange(0,101)

for i in range(len(voltage_range)):
    P=PEC_cell.PEC_Cell()
    P.Bandgap_message=False
    P.Eg_completed=False
    P.message_spectrum_loaded=False
    P.plot_curves=False
    P.set_voltage_max_ethylene=voltage_range[i]
    P.set_max_ethylene=100
    P.Calculate_setup(scenario='VariableCO2RRcatalyst')
    j_range[i]=np.max(P.Sparse_Coverage_J_matrix)
    print(r'At $\eta_{\rm CO2RR}$',voltage_range[i],'j_s=',j_range[i])
    print(P.Solar_cell_combination)

plt.figure()
plt.plot(voltage_range,j_range)
plt.xlabel(r'$\eta_{\rm CO_2RR}$ (V vs RHE)',fontsize=14)
plt.ylabel(r'Maximum $j_{\rm s}$ (mA/cm$^2$)',fontsize=14)




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
levels=[10,20,30,40,60,70]
clabel=r'C$_2$H$_4$ production rate ($\mu$mol/h/cm$^2$)'
figcontour, ax = plt.subplots()
  
    
CS1 = ax.contourf(X, Y, umol_production,levels=np.arange(-20,100,0.01),cmap=cmap)
CS2 = ax.contour(X, Y, umol_production,levels=levels,colors='w',linestyle='dashed')

plt.plot(voltage_range[np.array([0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,13])],np.array([68, 69, 71, 73, 75, 77, 79, 83, 85, 89, 92, 95.5, 97,100]),'--w')
#Cannot get a dashed contour plot for some reason                                        
                              
                                        
# CS3 = ax.contour(X, Y, umol_production,levels=[50],colors='w')
ax.clabel(CS2, fontsize=12,manual=[(-0.1,10),(-0.1,25),(-0.1,37),(-0.1,52),(-0.1,78),(-0.1,91.9)])
plt.plot(voltage_range,voltage_range*0+0.1,'w')

# cbar=figcontour.colorbar(CS1,ticks=np.arange(0,80,10))    
# cbar.ax.set_ylabel(clabel)

plt.xlabel(r'$\eta_{\rm CO_2RR}$ (V vs RHE)',fontsize=14)

plt.ylabel(r'$FE_{\rm C_2H_4}$ (%)',fontsize=14)

plt.plot(-0.61,61.2,'ok',markersize=9)
plt.plot(-0.53,58,'*k')

plt.ylim([0,99])
plt.yticks([0,20,40,60,80,100])

plt.plot(-1.55,60.7,'Xk',markersize=9)
plt.plot(-0.68,57,'^k',markersize=9)
plt.plot(-0.70,59.5,'Pk',markersize=9)
plt.plot(-1.05,51.5,'sk',markersize=9)
plt.plot(-1.05,38.7,'dk',markersize=9)
plt.plot(-1.00,33.2,'pk',markersize=9)

plt.show()
plt.figure()

boxtext='STE Efficiency: '
unit='%'
levels=np.arange(0,30,5)
cmap='Greens_r'
clabel='STE efficiency (%)'
figcontour, ax = plt.subplots()
    
CS1 = ax.contourf(X, Y, STE_production,levels=np.arange(-10,40,0.01),cmap=cmap)
CS2 = ax.contour(X, Y, STE_production,levels=levels,colors='w')
ax.clabel(CS2, fontsize=12,manual=[(-0.1,20),(-0.1,38),(-0.1,55),(-0.1,75),(-0.1,90)])
# cbar=figcontour.colorbar(CS1,ticks=levels)    
# cbar.ax.set_ylabel(clabel)

plt.xlabel(r'$\eta_{\rm CO_2RR}$ (V vs RHE)',fontsize=14)

plt.ylabel(r'$FE_{\rm C_2H_4}$ (%)',fontsize=14)

plt.plot(-0.61,61.2,'ok',markersize=9)
plt.plot(-0.53,58,'*k')

plt.plot(-1.55,60.7,'Xk',markersize=9)
plt.plot(-0.68,57,'^k',markersize=9)
plt.plot(-0.70,59.5,'Pk',markersize=9)
plt.plot(-1.05,51.5,'sk',markersize=9)
plt.plot(-1.05,38.7,'dk',markersize=9)
plt.plot(-1.00,33.2,'pk',markersize=9)

plt.show()

print('Isocurves')
