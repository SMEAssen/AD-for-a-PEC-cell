# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:08:36 2025

This script calculates the concentrated sunlight approach for PEC cells that reduce CO2 into C2H4 optimal at AM1.5G using a tandem solar cell. This generates Figure S2.  
Here, we use data from the O(II)D-Cu and NiFeOx catalysts, with R_solution=5 and FF solar cells=0.85

References
Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.
McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F


@author: assensme
"""

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'..\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell

B=PEC_cell.PEC_Cell()
B.message_spectrum_loaded=False
B.message_concentration=False
B.Calculate_setup(scenario='Scenario B: solar concentration')
plt.show()

print('Scenario B Solar concentration')