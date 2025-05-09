# -*- coding: utf-8 -*-
"""
Created on Wed May  7 11:47:31 2025

This script calculates the sparse coverage approach for PEC cells that reduce CO2 into C2H4 optimal at AM1.5G using a tandem solar cell and their response to surface degradation of 0-35%. This generates Figures S5 and S6
Here, we use data from the O(II)D-Cu, Cu-Ag and NiFeOx catalysts, with R_solution=5.

References
Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.
McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
Jing Gao  Dan Ren Xueyi Guo Shaik Mohammed Zakeeruddin and  Michael Grätzel  Faraday Discuss., 2019,215, 282-296

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