# -*- coding: utf-8 -*-
"""
Created on Sun May  4 14:23:58 2025

This script calculates the sparse coverage approach for PEC cells that reduce CO2 into C2H4 optimal at AM1.5G using a tandem solar cell and their response to variations in the sunlight. This generates Figure S9  
Here, we use data from the O(II)D-Cu, Cu-Ag and NiFeOx catalysts, with R_solution=5. The FF of the solar cells is 0.85, 0.75 and 0.65

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
A.message_spectrum_loaded=False
A.Calculate_setup(scenario='Scenario A: sparse coverage')
plt.show()

FF075=PEC_cell.PEC_Cell()
FF075.message_spectrum_loaded=False
FF075.Calculate_setup(scenario='FF 0.75 sparse coverage')
plt.show()

FF065=PEC_cell.PEC_Cell()
FF065.message_spectrum_loaded=False
FF065.Calculate_setup(scenario='FF 0.65 sparse coverage')
plt.show()

print('Scenario A sparse coverage')
