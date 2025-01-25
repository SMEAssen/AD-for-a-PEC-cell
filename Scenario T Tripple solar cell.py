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

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


#Tripple Solar cell
T=PEC_cell.PEC_Cell() 
T.plot_intersection=False
T.message_spectrum_loaded=False
T.Calculate_setup(scenario='Scenario T: sparse coverage with a triple junction')





 
# Creating dataset
z = np.arange(0.35,1.0,0.01)
y = np.arange(1.6,2.0,0.01)
x = np.arange(0.9,1.5,0.01)

X,Y,Z=np.meshgrid(x,y,z)
 
# Creating figure
fig = plt.figure(figsize = (16, 9))
ax = plt.axes(projection ="3d")
   
# Add x, y gridlines 
ax.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.3, 
        alpha = 0.2) 
 
 
# Creating color map
my_cmap = plt.get_cmap('jet')
 
# Creating plot
x_slice=[18,-1]



y_slice=[0]

z_slice=[34,0]

print(T.Sparse_Coverage_Efficiency_matrix[x_slice,y_slice,z_slice])

# ax.set_xlabel('Bandgap middle SC', fontweight ='bold') 
# ax.set_ylabel('Bandgap top SC', fontweight ='bold') 
# ax.set_zlabel('Bandgap bottom SC', fontweight ='bold')
# fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5,label='STE efficiency')
 
# show plot


fig = plt.figure(figsize = (16, 9))
ax = plt.axes(projection ="3d")
E=np.zeros(np.shape(T.Sparse_Coverage_Efficiency_matrix))

E[x_slice,:,:]=T.Sparse_Coverage_Efficiency_matrix[x_slice,:,:]+1
E[:,y_slice,:]=T.Sparse_Coverage_Efficiency_matrix[:,y_slice,:]+1
E[:,:,z_slice]=T.Sparse_Coverage_Efficiency_matrix[:,:,z_slice]+1


sctt = ax.scatter3D(X[E!=0], Y[E!=0], Z[E!=0],
                    c = E[E!=0]-1, 
                    cmap = my_cmap, vmax=16,vmin=0,
                    marker ='o')

ax.scatter3D(x[x_slice],y[y_slice],z[z_slice],'k') 
ax.set_xlabel('Bandgap middle SC', fontweight ='bold') 
ax.set_ylabel('Bandgap top SC', fontweight ='bold') 
ax.set_zlabel('Bandgap bottom SC', fontweight ='bold')
cbar=fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5,label='STE efficiency')




 

 
# Creating plot
x_slice=18
y_slice=29
z_slice=34


print(T.Sparse_Coverage_Efficiency_matrix[x_slice,y_slice,z_slice])

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
ax.set_xlabel('Bandgap middle SC', fontweight ='bold') 
ax.set_ylabel('Bandgap top SC', fontweight ='bold') 
ax.set_zlabel('Bandgap bottom SC', fontweight ='bold')
cbar=fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5,label='STE efficiency')



fig = plt.figure(figsize = (16, 9))
ax = plt.axes(projection ="3d")
   
# Add x, y gridlines 
ax.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.3, 
        alpha = 0.2) 
 
 
# Creating color map
my_cmap = plt.get_cmap('jet')
 
# Creating plot

sctt = ax.scatter3D(X, Y, Z,
                    c = T.Sparse_Coverage_Efficiency_matrix, 
                    cmap = my_cmap, 
                    marker ='o')
ax.scatter3D(1.15,1.75,0.55,'k')  
plt.title("simple 3D scatter plot")
ax.set_xlabel('X-axis', fontweight ='bold') 
ax.set_ylabel('Y-axis', fontweight ='bold') 
ax.set_zlabel('Z-axis', fontweight ='bold')
fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5)
 
# show plot
plt.show()

print('Tripple')
