# -*- coding: utf-8 -*-
"""
Created on Wed May  7 11:47:31 2025

@author: assensme
"""

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell
import numpy as np

OIID_Cu_Deg=PEC_cell.PEC_Cell()
OIID_Cu_Deg.message_spectrum_loaded=False
OIID_Cu_Deg.Calculate_setup(scenario='Scenario A: sparse coverage')



OIID_Cu_Deg.Degradation_catalyst(OIID_Cu_Deg.Solar_cell_combination[0],OIID_Cu_Deg.Solar_cell_combination[1],Degradation_range=np.arange(0,36,1))
plt.show()


Cu_Ag_Deg=PEC_cell.PEC_Cell()
Cu_Ag_Deg.message_spectrum_loaded=False
Cu_Ag_Deg.Calculate_setup(scenario='CuAg-NiFeOx_high_pH')


Cu_Ag_Deg.Degradation_catalyst(Cu_Ag_Deg.Solar_cell_combination[0],Cu_Ag_Deg.Solar_cell_combination[1],Degradation_range=np.arange(0,36,1))
plt.show()