# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:11:45 2025
This is a remake of Shu Hu et al. 

References
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
@author: assensme
"""

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

# #Check Hu H2 production
Hu=PEC_cell.PEC_Cell() 
Hu.Calculate_setup(scenario='Best_HER_Hu_2013')
plt.show()
print('Check Hu et al.')