# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:09:37 2025
This script calculates the PV-EC approach for reducing CO2 into C2H4 in an EC cell coupled to tandem solar cells at AM1.5G. This generates Figure S3.  
Here, we use data from the O(II)D-Cu and NiFeOx catalysts, with R_solution=5 and FF solar cells=0.85
First, the solar cells are calculated, which are coupled to external EC cells. In this configuration, 3 tandem solar cells could drive 2 EC cells.  

References
Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.
McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

@author: assensme
"""

import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'..\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell


C=PEC_cell.PEC_Cell()
C.message_spectrum_loaded=False
C.Bandgap_message=False
C.Calculate_setup(scenario='Scenario C: PV-EC')
plt.show()

print('Scenario C: PV-EC')
