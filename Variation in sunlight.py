# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:07:14 2025

This script calculates the sparse coverage approach for PEC cells that reduce CO2 into C2H4 optimal at AM1.5G using a tandem solar cell and their response to variations in the sunlight. This generates Figure 4 in the main text.  
Here, we use data from the O(II)D-Cu, Cu-Ag and NiFeOx catalysts, with R_solution=5 and FF solar cells=0.85

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


A=PEC_cell.PEC_Cell()
A.plot_curves=False
A.Bandgap_message=False
A.Eg_completed=False
A.message_spectrum_loaded=False
A.Calculate_setup(scenario='Scenario A: sparse coverage')
A.Vary_Solar_Intensity(A.Solar_cell_combination[0],A.Solar_cell_combination[1])
plt.show()

CuAg=PEC_cell.PEC_Cell()
CuAg.plot_curves=False
CuAg.Bandgap_message=False
CuAg.Eg_completed=False
CuAg.message_spectrum_loaded=False
CuAg.Calculate_setup(scenario='CuAg-NiFeOx_high_pH')
CuAg.Vary_Solar_Intensity(CuAg.Solar_cell_combination[0],CuAg.Solar_cell_combination[1])
plt.show()

print('Variation')