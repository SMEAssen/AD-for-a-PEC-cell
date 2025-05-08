# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:11:10 2025
This script calculates the sparse coverage approach for PEC cells that reduce CO2 into C2H4 optimal at AM1.5G using a triple solar cell. This generates Figure 3 in the main text.  
Here, we use data from the O(II)D-Cu and NiFeOx catalysts, with R_solution=5 and FF solar cells=0.85

References
Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.
McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

@author: assensme
"""

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'...\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell
import numpy as np



#Tripple Solar cell
T=PEC_cell.PEC_Cell() 
T.plot_intersection=False
T.plot_curves=False
T.message_spectrum_loaded=False
T.Calculate_setup(scenario='Scenario T: sparse coverage with a triple junction')


print('Max Efficiency',np.max(T.Sparse_Coverage_Efficiency_matrix), 'with solar cells', T.Solar_cell_combination)
 
# Creating dataset
z = np.arange(0.35,1.0,0.01)
y = np.arange(1.6,2.0,0.01)
x = np.arange(0.9,1.5,0.01)

X,Y,Z=np.meshgrid(x,y,z)
 
 

 
# Creating plot
x_slice=18
y_slice=29
z_slice=34





fig = plt.figure(figsize = (16, 9))
ax = plt.axes(projection ="3d")
E=np.zeros(np.shape(T.Sparse_Coverage_Efficiency_matrix))

E[x_slice,:,:]=T.Sparse_Coverage_Efficiency_matrix[x_slice,:,:]+1
E[:,y_slice,:]=T.Sparse_Coverage_Efficiency_matrix[:,y_slice,:]+1
E[:,:,z_slice]=T.Sparse_Coverage_Efficiency_matrix[:,:,z_slice]+1

# Creating color map
my_cmap = plt.get_cmap('jet')
sctt = ax.scatter3D(X[E!=0], Y[E!=0], Z[E!=0],
                    c = E[E!=0]-1, 
                    cmap = my_cmap, vmax=16,vmin=0,
                    marker ='o')

ax.scatter3D(x[x_slice],y[y_slice],z[z_slice],'k') 
ax.set_xlabel(r'$E_{\rm g2}$' ,fontsize=16) 
ax.set_ylabel(r'$E_{\rm g1}$' ,fontsize=16) 
ax.set_zlabel(r'$E_{\rm g3}$',fontsize=16)
cbar=fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5)
cbar.ax.set_ylabel(r'$\eta_{\rm STE}$ (%)',fontsize=16)




 
# show plot
plt.show()

print('Scenario T: Triple Solar cell, sparse coverage')


X,Y=np.meshgrid(x,y)

plt.figure()
figcontour, ax = plt.subplots()


Z=E[:,:,z_slice]=T.Sparse_Coverage_Efficiency_matrix[:,:,z_slice]+1
levels=np.arange(0,18,2)

cmap='jet'
CS1 = ax.contourf(X, Y, Z,levels=np.arange(0,16.01,0.01),cmap=cmap)  

cbar=figcontour.colorbar(CS1,ticks=levels)    
clabel=r'$\eta_{\rm STE}$ (%)'

cbar.ax.set_ylabel(clabel,fontsize=14)
xmax=0
ymax=0
if xmax==0 and ymax==0:
    xmax=np.where(Z==np.max(Z))[1][0]
    ymax=np.where(Z==np.max(Z))[0][0]
    
plt.plot(x[xmax],y[ymax],'ok')
ax.set_xlabel(r'$\it{E}$$_{\rm g2}$ (eV)',fontsize=14)
ax.set_ylabel(r'$\it{E}$$_{\rm g1}$ (eV)',fontsize=14)
props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=1)

plt.show()

X,Y=np.meshgrid(z,x)

plt.figure()
figcontour, ax = plt.subplots()


Z=T.Sparse_Coverage_Efficiency_matrix[x_slice,:,:]+1
levels=np.arange(0,18,2)

cmap='jet'
CS1 = ax.contourf(X, Y, Z,levels=np.arange(0,16.01,0.01),cmap=cmap)  
cbar=figcontour.colorbar(CS1,ticks=levels)    
clabel=r'$\eta_{\rm STE}$ (%)'

cbar.ax.set_ylabel(clabel,fontsize=14)
xmax=0
ymax=0
if xmax==0 and ymax==0:
    xmax=np.where(Z==np.max(Z))[1][0]
    ymax=np.where(Z==np.max(Z))[0][0]
    
plt.plot(z[xmax],x[ymax],'ok')
ax.set_xlabel(r'$\it{E}$$_{\rm g3}$ (eV)',fontsize=14)
ax.set_ylabel(r'$\it{E}$$_{\rm g2}$ (eV)',fontsize=14)
props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=1)
  
plt.show()


X,Y=np.meshgrid(z,y)

plt.figure()
figcontour, ax = plt.subplots()


Z=T.Sparse_Coverage_Efficiency_matrix[:,y_slice,:]+1
levels=np.arange(0,18,2)

cmap='jet'
CS1 = ax.contourf(X, Y, Z,levels=np.arange(0,16.01,0.01),cmap=cmap)  
cbar=figcontour.colorbar(CS1,ticks=levels)    
clabel=r'$\eta_{\rm STE}$ (%)'

cbar.ax.set_ylabel(clabel,fontsize=14)
xmax=0
ymax=0
if xmax==0 and ymax==0:
    xmax=np.where(Z==np.max(Z))[1][0]
    ymax=np.where(Z==np.max(Z))[0][0]
    
plt.plot(z[xmax],y[ymax],'ok')
ax.set_xlabel(r'$\it{E}$$_{\rm g3}$ (eV)',fontsize=14)
ax.set_ylabel(r'$\it{E}$$_{\rm g1}$ (eV)',fontsize=14)
props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=1)

plt.show()
print('Scenario T: Triple Solar cell, sparse coverage')
