# -*- coding: utf-8 -*-
"""
Created on Sun May  4 14:23:58 2025

@author: assensme
"""
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(r'C:\Users\assensme\Documents\Programming\PEC\PEC_cell.py')))
import matplotlib.pyplot as plt
import PEC_cell

A=PEC_cell.PEC_Cell()
A.message_spectrum_loaded=False
A.Calculate_setup(scenario='Scenario A: sparse coverage')
plt.show()

B=PEC_cell.PEC_Cell()
B.message_spectrum_loaded=False
B.Calculate_setup(scenario='FF 0.75 sparse coverage')
plt.show()

C=PEC_cell.PEC_Cell()
C.message_spectrum_loaded=False
C.Calculate_setup(scenario='FF 0.65 sparse coverage')
plt.show()

print('Scenario A sparse coverage')
