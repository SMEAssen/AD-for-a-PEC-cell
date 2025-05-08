# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:28:41 2023

This script belongs to the paper "Axiomatic Design for Analysis and Noniterative Engineering of Photoelectrochemical Cells for CO2 Conversion at High Photon to Product Yield"

Four types of concentration are distinguished: sparse coverage, solar concentration, PV-EC and no concentration, each with their own module. 
The reduction reaction is assumed to be CO2RR towards ethylene, but the model can also be used for hydrogen evolution reaction (HER) calculations. 
Different CO2RR catalysts for ethylene production are parameterized in the required voltage, current density and percentage ethylene, similar for OER catalysts.
Additionally, the fluid resistance (standard: 5Ω) and FF (standard 0.85) are specified. 

Using these parameters and the AM1.5G or AM1.5D spectrum23, V, j_s,FE_C2H4,ν_C2H4 and η_STE are calculated for different combinations of tandem solar cells. 
The high bandgap solar cell is assumed to be transparent for all wavelengths above its absorption threshold, allowing these photons to be absorbed by a medium bandgap solar cell. 
For each available current density, the voltages of the first and second solar cells are added together, resulting in a tandem solar performance graph similar to Figure 2b (main text). 
The intersection of the tandem solar cells graph with the catalyst graph gives j_s, used to calculate V, j_s,FE_C2H4 (based on interpolated data), ν_C2H4 and η_STE. 
For the OER catalyst, it is assumed j_s=j_OER, and J=j_s⋅1cm2. 

For the scenarios where the optimal solar cells are calculated with sparse coverage, the CO2RR catalyst is assumed to operate at the overpotential that generates the highest Faradaic efficiency. 
For every solar cell combination, a different j_s is calculated, which in turn is used to determine the concentration ratio A_s/A_(CO_2 RR) for the sparsely applied catalyst scenario. 
In the concentrated sunlight scenario, the solar concentration on the solar cells is gradually increased to achieve the ideal j_CO2RR=j_s, optimizing the selectivity towards ethylene. 
Additionally, R_transport is lowered to 1Ω to account for the larger surface area of the CO2RR catalyst. In the PV-EC scenario, all current densities are calculated separately. 
The optimal solar cells for electricity production are calculated and matched to an EC cell that operates at a fixed voltage. 
The scaling, both of the different area and number of PV and EC components, is determined afterwards. 
When varying sunlight and calculating the scenario without concentration, the Faradaic  efficiency, voltage and current responses of the CO2RR catalysts are interpolated from reported data.

In this file, 5 calculations are preselected.  
Using scenario A: sparse coverage gives Figure 3. The combination of solar cells yielding the higest ν_C2H4 is used to calculate variations in solar light gives Figure 4a. 
Figure 4b is made by calculating the bandgaps with the largest production rate in the CuAg-NiFeOx_high_pH scenario and varying the solar irradiation. 
Figure S2 is made with ‘Scenario B: solar concentration’, Figure S3 with‘Scenario C: PV-EV’ and Figure S6 with ‘Scenario F: no concentration’.    


References

Calculations:
C. H. Henry, Journal of Applied Physics 51, 4494 (1980); https://doi.org/10.1063/1.328272
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
Rühle, S. Tabulated Values of the Shockley–Queisser Limit for Single Junction Solar Cells. Solar Energy 2016, 130, 139–147; https://doi.org/10.1016/j.solener.2016.02.015.


CO2RR to ethylene catalysts:
Dan Ren ACS Sustainable Chemistry & Engineering,5(10):9191–9199 (2017); https://doi.org/10.1021/acssuschemeng.7b02110
Cao-Thang Dinh et al. ,CO2 electroreduction to ethylene via hydroxide-mediated copper catalysis at an abrupt interface.Science360,783-787(2018).DOI:10.1126/science.aas9100
Ren D, Gao J, Zakeeruddin SM, Grätzel M. New Insights into the Interface of Electrochemical Flow Cells for Carbon Dioxide Reduction to Ethylene. J Phys Chem Lett. 2021 Aug 12;12(31):7583-7589. doi: 10.1021/acs.jpclett.1c02043. Epub 2021 Aug 4. PMID: 34347495.
Jing Gao, Hong Zhang, Xueyi Guo, Jingshan Luo, Shaik M. Zakeeruddin, Dan Ren*, and Michael Grätzel*J. Am. Chem. Soc. 2019, 141, 47, 18704–18714
Jing Gao  Dan Ren Xueyi Guo Shaik Mohammed Zakeeruddin and  Michael Grätzel  Faraday Discuss., 2019,215, 282-296
Dongxing Tan, Bari Wulan, Xueying Cao, Jintao Zhang, Nano Energy,Volume 89, Part B, 2021, 106460, 2211-2855 https://doi.org/10.1016/j.nanoen.2021.106460.
Choi, C., Kwon, S., Cheng, T. et al. Highly active and stable stepped Cu surface for enhanced electrochemical CO2 reduction to C2H4. Nat Catal 3, 804–812 (2020). https://doi.org/10.1038/s41929-020-00504-x
Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.

HER catalyst:
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

OER catalyst:
Shu H et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
Joya & de Groot ACS Catal. 2016, 6, 3, 1768–1771
McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  



Solar Spectrum: https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html
Intersect:https://github.com/sukhbinder/intersection/blob/master/intersect/intersect.py

@author: assensme
"""


import numpy as np
import matplotlib.pyplot as plt
import copy

import time
from scipy import interpolate
from scipy.special import lambertw


class intersect():
    #Source: https://github.com/sukhbinder/intersection/blob/master/intersect/intersect.py
    def _rect_inter_inner(x1, x2):
        n1 = x1.shape[0]-1
        n2 = x2.shape[0]-1
        X1 = np.c_[x1[:-1], x1[1:]]
        X2 = np.c_[x2[:-1], x2[1:]]
        S1 = np.tile(X1.min(axis=1), (n2, 1)).T
        S2 = np.tile(X2.max(axis=1), (n1, 1))
        S3 = np.tile(X1.max(axis=1), (n2, 1)).T
        S4 = np.tile(X2.min(axis=1), (n1, 1))
        return S1, S2, S3, S4
    
    
    def _rectangle_intersection_(x1, y1, x2, y2):
        S1, S2, S3, S4 = intersect._rect_inter_inner(x1, x2)
        S5, S6, S7, S8 = intersect._rect_inter_inner(y1, y2)
    
        C1 = np.less_equal(S1, S2)
        C2 = np.greater_equal(S3, S4)
        C3 = np.less_equal(S5, S6)
        C4 = np.greater_equal(S7, S8)
    
        ii, jj = np.nonzero(C1 & C2 & C3 & C4)
        return ii, jj
    
    
    def intersection(x1, y1, x2, y2):
        x1 = np.asarray(x1)
        x2 = np.asarray(x2)
        y1 = np.asarray(y1)
        y2 = np.asarray(y2)
    
        ii, jj = intersect._rectangle_intersection_(x1, y1, x2, y2)
        n = len(ii)
    
        dxy1 = np.diff(np.c_[x1, y1], axis=0)
        dxy2 = np.diff(np.c_[x2, y2], axis=0)
    
        T = np.zeros((4, n))
        AA = np.zeros((4, 4, n))
        AA[0:2, 2, :] = -1
        AA[2:4, 3, :] = -1
        AA[0::2, 0, :] = dxy1[ii, :].T
        AA[1::2, 1, :] = dxy2[jj, :].T
    
        BB = np.zeros((4, n))
        BB[0, :] = -x1[ii].ravel()
        BB[1, :] = -x2[jj].ravel()
        BB[2, :] = -y1[ii].ravel()
        BB[3, :] = -y2[jj].ravel()
    
        for i in range(n):
            try:
                T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
            except:
                T[:, i] = np.Inf
    
        in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (
            T[0, :] <= 1) & (T[1, :] <= 1)
    
        xy0 = T[2:, in_range]
        xy0 = xy0.T
        return xy0[:, 0], xy0[:, 1] 
    


class PEC_Cell():
    
    def __init__(self):
        self.time_init=time.time()
        self.Solar_spectrum_location   = r"C:\Users\assensme\Documents\PhD\Modelling\Solar spectrum\AM1.5.csv" #Location of data with AM1.5 G 1307 data, taken from nrel site. 
     
        self.V_OER=1.23
        self.V_CO2R=0.08
        self.V_HER_potential=0.00
        
        self.Temperature=298.15 #K
        
        self.nmeV=1239.84207   # nm to eV 
        
        self.h=6.62607015e-34 #m^2 kg/s
        self.c=299792458 #m/s
        self.A=6.0221409e+23 #Avogadro
        self.q=1.602176634e-19 #J/eV
        self.F=self.A*self.q
        self.kb=1.38064852e-23 #m^2kg/s^2/K
        self.R=8.3145 #J/mol/K
        
        self.kTeV           = 0.02585202874091 * self.Temperature / 300.0       # kT in eV
        self.kTJ=self.kTeV * self.q 
        
        
        self.Second_solar_cell=False
        self.Three_solar_cells=False
        
        
        self.PEC_mode='sparse coverage' #Scenario A
        self.Fixed_ethylene_output=True
        
        self.Reduction_mode='CO2 reduction'
        
        self.HER_catalyst=''   
        self.CO2RR_catalyst='OIID-Cu'
        self.OER_catalyst='NiFeOx_pH14'
        
        
        self.Fluid_resistance=5 #Ohm/cm^2
        self.Voltage_accuracy=0.001 #V
        self.Current_accuracy=0.01 # mA
        
        
        self.Spectrum_Loaded=False
        self.Mirror_factor=1
        self.Solar_Concentration=1
        self.FF_goal=0.85
        self.Co2RcatalystConcentrator=1
        
        self.bandgap_range=np.array([1.6,2.4,1.0,1.9])
        self.grid_step=np.array([0.05,0.05])
        
        self.plot_spectrum=False
        self.plot_intersection=False
        
        self.Bandgap_message=True
        self.message_spectrum_loaded=True
        self.Eg_completed=True
        self.plot_curves=True
        self.message_concentration=True
        self.message_variation=True
        
    def Calculate_setup(self,scenario='OD(II)-Cu_NiFeOx_high_pH_sparse_coverage'):
        #Start up a scenario and catalyst set-up and start the calculations 
        self.Set_scenario(scenario)
        self.Set_up_catalyst()
        
            
        if self.PEC_mode=='Solar cells':
            self.Calculate_solar_cells()
        
        elif self.PEC_mode=='1:1':
            self.Calculate_1on1()
        
        elif self.PEC_mode=='sparse coverage':
            self.Calculate_sparse_coverage()
            
        elif self.PEC_mode=='solar concentration':
            self.Calculate_solar_light_concentration()
            
        elif self.PEC_mode=='PV-EC':
            self.Calculate_solar_cells()
            self.Calculate_converter_PV_EC()


        
    def Set_scenario(self,scenario):
        #Different scenarios
        
        if scenario=='Best_HER_Hu_2013':
            #Remake of Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
            self.PEC_mode='1:1' #No concentration
            self.Reduction_mode='HER'
            self.HER_catalyst='Platinum'
            self.OER_catalyst='RuO2 Neutral'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
            
        elif scenario=='Common_earth_HER_Hu_2013':
            #Remake of Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
            self.PEC_mode='1:1'
            self.Reduction_mode='HER'
            self.HER_catalyst='Common Earth, Ni-MO'
            self.OER_catalyst='NiFeOx Neutral'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
            
        elif scenario=='Scenario F: no concentration':
            #Failed PEC cell for CO2RR, without any concentration
            self.PEC_mode='1:1' #No concentration
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='Scenario A: sparse coverage':
            # print('A')
  
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm
            self.FF_goal=0.85
            
            
        elif scenario=='FF 0.75 sparse coverage':
            # print('A')
  
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm
            self.FF_goal=0.75
                
        elif scenario=='FF 0.65 sparse coverage':
            # print('A')
  
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm
            self.FF_goal=0.65
            
        elif scenario=='Scenario B: solar concentration':
         
            self.PEC_mode='solar concentration'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=1 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='Scenario C: PV-EC':
         
            self.PEC_mode='PV-EC'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='CuAg-NiFeOx_high_pH':
            #Different catalyst, used for the variance in solar light graph
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='CuAg'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2  
            self.FF_goal=0.85
            
            
            
        elif scenario=='Scenario T: sparse coverage with a triple junction':
            #Scenario A with a triple junction
            self.Three_solar_cells=True
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85    
            
        elif scenario=='CO2RR_Greatzel2022_NiFeOx_high_pH':
            #Similar to Scenario A, with a different (better) CO2RR catalyst, not plublished yet
            self.PEC_mode='sparse coverage' 
            self.Reduction_mode='HER'
            self.CO2RR_catalyst='Greatzel2022'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='OD(I)-Cu_NiFeOx_high_pH_sparse_coverage':
            #Different catalyst
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='Solar cells & Electrolyzers OD(II)-Cu_NiFeOx_high_pH':
            self.PEC_mode='Solar cells'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='OIID-Cu'
            self.OER_catalyst='NiFeOx_pH14'
            self.FF_goal=0.85
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='Variable'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        elif scenario=='VariableCO2RRcatalyst':
            self.PEC_mode='sparse coverage'
            self.Reduction_mode='CO2 reduction'
            self.CO2RR_catalyst='Variable'
            self.OER_catalyst='NiFeOx_pH14'
            self.Fluid_resistance=5 #Ohm/cm^2
            self.FF_goal=0.85
            
        #The resolution of the graphs is determined here    
        if self.PEC_mode=='Solar cells' or self.PEC_mode=='PV-EC':
            self.bandgap_range=np.array([1.35,2.4,0.5,1.3])
            self.grid_step=np.array([0.01,0.01])
            
        elif self.Reduction_mode=='HER': 
            self.bandgap_range=np.array([1.35,2.45,0.65,1.325])
            self.grid_step=np.array([0.05,0.025])
            

        elif self.Reduction_mode=='CO2 reduction':
            if self.Three_solar_cells:
                
                self.bandgap_range=np.array([1.6,2.0,0.9,1.5,0.35,1.0])
                self.grid_step=np.array([0.01,0.01,0.01])   
            else:
                self.bandgap_range=np.array([1.6,2.6,1.0,2.1])
                self.grid_step=np.array([0.05,0.05])   
        

    def Set_up_catalyst(self):
        #Gathers data for OER and HER/CO2RR catalyst
        
        
        self.i_OER,self.eta_OER,self.A_OER=self.Data_OER()  
        if self.Reduction_mode=='CO2 reduction':
            if self.CO2RR_catalyst=='Variable':
                self.V_max_ethylene=self.set_voltage_max_ethylene 
                self.J_max_ethylene=250
                self.Max_ethylene_percentage=self.set_max_ethylene
            else:
                
                Ethylene_percentage,Current_density_CO2RR,Voltage_CO2RR=self.Data_CO2RR_catalyst()
                self.eta_CO2RR=np.arange(max(Voltage_CO2RR),min(Voltage_CO2RR)+self.Voltage_accuracy,-self.Voltage_accuracy)
                CD_extended=interpolate.interp1d(Voltage_CO2RR,Current_density_CO2RR)
                self.Current_density_CO2RR_extended=CD_extended(self.eta_CO2RR)
                Ethylene_extended_function=interpolate.interp1d(Voltage_CO2RR,Ethylene_percentage,'cubic')
                
                self.Ethylene_percentage_extended=Ethylene_extended_function(self.eta_CO2RR)
                #self.Ethylene_percentage_extended[np.where(self.Ethylene_percentage_extended<min(Ethylene_percentage))]=min(Ethylene_percentage)
                location_max_ethylene=np.where(Ethylene_percentage==np.max(Ethylene_percentage))[0][0]
                self.V_max_ethylene=Voltage_CO2RR[location_max_ethylene]       
                self.J_max_ethylene=Current_density_CO2RR[location_max_ethylene]
                self.Max_ethylene_percentage=np.max(Ethylene_percentage)
                # print(self.Max_ethylene_percentage,self.J_max_ethylene,self.V_max_ethylene,self.CO2RR_catalyst)
        elif self.Reduction_mode=='HER':
            self.i_HER,self.eta_HER,self.A_HER=self.Data_HER()

    def Data_CO2RR_catalyst(self):
        if self.CO2RR_catalyst=='DanRen2017':
            #ACS Sustainable Chem. Eng. 2017, 5, 10, 9191–9199
            
            Ethylene_percentage=[4,14.5,26.8,30,24,20.5]
            Current_density=np.array([70,47,36,25.5,19,11])
            Voltage=np.array([-1.1,-1.05,-1,-0.95,-0.9,-0.85])
            
        elif self.CO2RR_catalyst=='Dinh2018':
            #Cao-Thang Dinh et al. ,CO2 electroreduction to ethylene via hydroxide-mediated copper catalysis at an abrupt interface.Science360,783-787(2018).DOI:10.1126/science.aas9100
           
            Ethylene_percentage=[65,70,67,58,50,40,0]
            Current_density=[150,100,60,35,25,15,10]
            Voltage=[-0.58,-0.57,-0.53,-0.47,-0.41,-0.32,-0.25]
            
            
        elif self.CO2RR_catalyst=='Greatzel2022':
            # Sun2Chem presentation WP4 2022
            Ethylene_percentage=[19,52,63,75,60,30]
            Current_density=[10,30,70,100,150,200]
            Voltage=[-0.35,-0.52,-0.56,-0.58,-0.615,-0.64]
         
        
        elif self.CO2RR_catalyst=='Cu200nm':
            #Ren D, Gao J, Zakeeruddin SM, Grätzel M. New Insights into the Interface of Electrochemical Flow Cells for Carbon Dioxide Reduction to Ethylene. J Phys Chem Lett. 2021 Aug 12;12(31):7583-7589. doi: 10.1021/acs.jpclett.1c02043. Epub 2021 Aug 4. PMID: 34347495.
            
            Ethylene_percentage=[0.03,0.43,2.06,11.05,22.3,43.25,51.21,59.49,34.56]
            Current_density=[18.92,24.31,33.38,32.15,36.88,94.06,181.37,282.73,390]
            Voltage=np.arange(-0.35,-0.80,-0.05)
       
            
        elif self.CO2RR_catalyst=='OIID-Cu':
            #Ren D, Gao J, Zakeeruddin SM, Grätzel M. New Insights into the Interface of Electrochemical Flow Cells for Carbon Dioxide Reduction to Ethylene. J Phys Chem Lett. 2021 Aug 12;12(31):7583-7589. doi: 10.1021/acs.jpclett.1c02043. Epub 2021 Aug 4. PMID: 34347495.
            
            Ethylene_percentage=[17,17.25,43,57.23,58.73,61.18,59.79,57.28,52.13]
            Current_density=[20,50,100,150,200,250,300,350,400]
            Voltage=[-0.40,-0.44,-0.5375,-0.575,-0.60,-0.615,-0.627,-0.6375,-0.65]
        
        elif self.CO2RR_catalyst=='OID-Cu':
            #Ren D, Gao J, Zakeeruddin SM, Grätzel M. New Insights into the Interface of Electrochemical Flow Cells for Carbon Dioxide Reduction to Ethylene. J Phys Chem Lett. 2021 Aug 12;12(31):7583-7589. doi: 10.1021/acs.jpclett.1c02043. Epub 2021 Aug 4. PMID: 34347495.
            
            Ethylene_percentage=[0,14.02,42.75,49.02,50.74,55.94,53.25,53.09,52.63]
            Current_density=[20,50,100,150,200,250,300,350,400]
            Voltage=[-0.60,-0.61,-0.66,-0.68,-0.70,-0.72,-0.73,-0.74,-0.75]           
            
        elif self.CO2RR_catalyst=='Tan2021nearNeutral':
            #Dongxing Tan, Bari Wulan, Xueying Cao, Jintao Zhang, Nano Energy,Volume 89, Part B, 2021, 106460, 2211-2855 https://doi.org/10.1016/j.nanoen.2021.106460.

            Ethylene_percentage=[40,47,55,60,69.8,55,41,33,33]
            Current_density=[37.5,33,29.5,25.5,22.5,19,16.5,15,12]
            Voltage=np.arange(-1.3,-0.85,0.05)
            
        
        elif self.CO2RR_catalyst=='Tan2021KOH':
            #Dongxing Tan, Bari Wulan, Xueying Cao, Jintao Zhang, Nano Energy,Volume 89, Part B, 2021, 106460, 2211-2855 https://doi.org/10.1016/j.nanoen.2021.106460.


            Voltage=np.array([-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6])
            Current_density=np.array([170,155,140,130,120,112,105,102])
            Ethylene_percentage=np.array([38,48,52.5,59.5,66,76,62.5,48.5])
       
            
        elif self.CO2RR_catalyst=='CuAu':
            #Faraday Discuss., 2019,215, 282-296
            Voltage=np.array([-0.6,-0.7,-0.8,-0.9,-0.95,-1,-1.05,-1.1,-1.15])
            Current_density=np.array([2,4.53,7.7,17.73,24.32,31.24,42.87,48.82,63.98])
            Ethylene_percentage=np.array([0,1.15,4.8,31.7,34.1,36.17,38.7,34.97,15.85])
 
            
        elif self.CO2RR_catalyst=='CuAg':
            #J. Am. Chem. Soc. 2019, 141, 47, 18704–18714
            Voltage=np.array([-0.6,-0.7,-0.8,-0.9,-0.95,-1,-1.05,-1.1,-1.15,-1.2,-1.25,-1.3])
            Current_density=np.array([1.82,3.69,6.2,12.66,18.62,28.58,35.12,48.04,64.43,80.29,110,200])
            Ethylene_percentage=np.array([0.57,0.93,3.30,12.93,24.00,34.03,51.50,35.33,26.13,15.70,10,5])
    
        elif self.CO2RR_catalyst=='Copper':
            
            #J. Am. Chem. Soc. 2019, 141, 47, 18704–18714
            Voltage=np.array([-0.6,-0.7,-0.8,-0.9,-0.95,-1,-1.05,-1.1,-1.15,-1.2])
            Current_density=np.array([2.23,3.95,9.13,12.85,19,25.47,28.08,35.2,54.43,68.94])
            Ethylene_percentage=np.array([0,1.8,4.17,10.80,23.83,33.20,30.47,25.47,24.07,12.07])
          
        elif self.CO2RR_catalyst=='Choi2020':
            #Choi, C., Kwon, S., Cheng, T. et al. Highly active and stable stepped Cu surface for enhanced electrochemical CO2 reduction to C2H4. Nat Catal 3, 804–812 (2020). https://doi.org/10.1038/s41929-020-00504-x
            Voltage=np.array([-0.8,-0.95,-0.98,-1,-1.03,-1.07])
            Voltage_partial_CD=np.array([-0.78,-0.82,-0.9,-0.95,-0.98,-1.01,-1.06,-1.1])
            Partial_Current_density_original=np.array([0,0.5,1,4,15,23,25,28]) #Partial current density towards Ethyelne
            Ethylene_percentage=np.array([15,57,66,72,69,61])
            
            cd=interpolate.interp1d(Voltage_partial_CD,Partial_Current_density_original)
            Current_density=cd(Voltage)/Ethylene_percentage*100
        else:
            print('No CO2RR catalyst selected')
            
        
        return Ethylene_percentage,Current_density,Voltage
        
    def Data_HER(self):
        #Shu Hu, Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
        
        if self.HER_catalyst=='Platinum':
            i_HER=10 #mA/cm^2
            eta_HER=55e-3 #V
            A_HER=30e-3 #V/dec, Tafel slope
        
        elif self.HER_catalyst=='Common Earth, Ni-MO':
            i_HER=10 #mA/cm^2
            eta_HER=75e-3 #V
            A_HER=40e-3 #V/dec, Tafel slope
        else:
            print('No HER catalyst selected')
        return i_HER,eta_HER,A_HER
    
    def Data_OER(self):
        
    
        if self.OER_catalyst=='RuO2 Neutral':
            #Shu Hu, Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
            i_OER=10 #mA/cm^2
            eta_OER=240e-3 #V
            A_OER=37e-3 #V/dec, Tafel slope
        elif self.OER_catalyst=='NiFeOx Neutral':
            #Shu Hu, Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F
            i_OER=10 #mA/cm^2
            eta_OER=280e-3 #V
            A_OER=40e-3 #V/dec, Tafel slope
            
        elif self.OER_catalyst=='CuxO':
            #Joya & de Groot ACS Catal. 2016, 6, 3, 1768–1771
            i_OER=50 #mA/cm^2
            eta_OER=1.84-1.23 #V
            A_OER=44e-3 #V/dec, Tafel slope 
            
        elif self.OER_catalyst=='NiFeOx_pH14':
            #McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.
            i_OER=10#mA/cm^2
            eta_OER=327e-3#V
            A_OER=34e-3#V/dec, Tafel slope 
            
        elif self.OER_catalyst=='IrOx_pH14':
            #McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.
            i_OER=10#mA/cm^2
            eta_OER=325e-3#V
            A_OER=40e-3#V/dec, Tafel slope 
            
            
        elif self.OER_catalyst=='NiCoOx_pH14':
            #McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.
            i_OER=10#mA/cm^2
            eta_OER=345e-3#V
            A_OER=33e-3#V/dec, Tafel slope 
            
            
        
        else: 
            print('No OER catlayst selected')
        
        return i_OER,eta_OER,A_OER
      
 
    
    def Calculate_solar_cells(self):
        #This calculates the solar cell efficiency for conversion towards electricity. It is only used for calculations with an electrolyzer later on
        
        time_start=time.time()
        #Load the spectrum and set up the matrices
        self.Load_Spectrum_AM15()
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        if self.Three_solar_cells:
            Solar_cell_3_range=np.arange(self.bandgap_range[4],self.bandgap_range[5],self.grid_step[2])
        else:
            Solar_cell_3_range=[0]
        Efficiency_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        V_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        J_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        
        #Calculations for every solar cell
        count_x=0
        for Bandgap_1 in Solar_cell_1_range:
            count_y=0
            
            for Bandgap_2 in Solar_cell_2_range:
                
                count_z=0
                for Bandgap_3 in Solar_cell_3_range:
                
                    #The solar cell with the bigger bandgap should go on top
                    if Bandgap_2>=Bandgap_1 or Bandgap_3>Bandgap_1 or Bandgap_3>Bandgap_2:
                        continue
                    #The solar cells tandem efficiency is caluclated here
                    if self.Three_solar_cells:
                        
                        Voltage_SC1,J_SC1,FF,Voltage_SC2,J_SC2,FF,Voltage_SC3,J_SC3,FF=self.Triple_solar_cell(Bandgap_1,Bandgap_2,Bandgap_3)
                        V_tandem_solar_cell,J_tandem_solar_cell=self.find_J_V_triple_solar_cel(Voltage_SC1,J_SC1,Voltage_SC2,J_SC2,Voltage_SC3,J_SC3)
                        
                    else:
                        Voltage_SC1,J_SC1,FF,Voltage_SC2,J_SC2,FF=self.Double_solar_cell(Bandgap_1,Bandgap_2)
                        V,J=self.find_J_V_tandem_solar_cel(Voltage_SC1,J_SC1,Voltage_SC2,J_SC2)
                    #Store data
                    eta,V_matrix[count_x,count_y,count_z],J_matrix[count_x,count_y,count_z]=self.Efficiency_solar_cell(V,J)
                    Efficiency_matrix[count_x,count_y,count_z]=eta*100

                    if self.Bandgap_message:
                        print(Bandgap_1,Bandgap_2,Efficiency_matrix[count_x,count_y,count_z])
                    count_z+=1
                
                count_y+=1
            count_x+=1
            if self.Eg_completed:
                print('Time so far:' +'{0:.2f}'.format((time.time()-time_start)/60)+'Minutes')
        #Store data    
        self.Max_efficiency_Solar_cell_3_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[2][0]    
        self.Max_efficiency_Solar_cell_2_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[1][0]
        self.Max_efficiency_Solar_cell_1_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[0][0]
        self.Solar_cell_combination=[Solar_cell_1_range[self.Max_efficiency_Solar_cell_1_index],Solar_cell_2_range[self.Max_efficiency_Solar_cell_2_index],Solar_cell_3_range[self.Max_efficiency_Solar_cell_3_index]]
        
        self.Efficiency_matrix=Efficiency_matrix
        self.Solar_V_matrix=V_matrix
        self.Solar_J_matrix=J_matrix
        
        if not self.Three_solar_cells and self.plot_curves:
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,Efficiency_matrix[:,:,0],xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext='Efficiency: ',unit='%',levels=[0,10,15,20,25,30,35,40,45,50],detail='%1.1f',cmap='Greens_r',title='Solar cell tandem series efficiency',clabel='Solar cell Efficiency (%)',decimals=1)
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,V_matrix[:,:,0],xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext=r'$U_{\rm solar}$: ',title='',unit='V',levels=np.arange(np.min(V_matrix),np.max(V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',clabel=r'$U_{\rm solar}$: (V)',decimals=2)
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,J_matrix[:,:,0],xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext=r'$j_{\rm s}$',title='',unit=' mA/cm$^2$',cmap='autumn',levels=np.array([0,3,6,9,12,15,18,21,24,27,30,33,36]),clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
           
 
    def Calculate_converter_PV_EC(self,eta_membrane=0.1,converter_efficiency=1):
        #Using solar cells at the MPP, calculate the PV-EC efficiency using an EC module. We assume that the catalyst operates at maximum faradaic efficiency. 
        
        #Determine electrolyzer efficiency
        eta_OER=self.V_OER+self.eta_OER
        eta_CO2RR=self.V_max_ethylene
        electrolyzer_efficiency=(self.V_OER-self.V_CO2R)/(eta_OER-eta_CO2RR+eta_membrane)
        #print(electrolyzer_efficiency)
        
        Faradaic_efficiency=self.Max_ethylene_percentage/100
        
        #J on the electrolyzers
        j_op_converted=self.Efficiency_matrix/100*self.Solar_Power/(eta_OER-eta_CO2RR+eta_membrane)
        j_C2H4=j_op_converted*Faradaic_efficiency
        
        #calculate efficiencies
        self.product_solar_efficiency=j_C2H4*1.15/self.Solar_Power
        self.umol_cm2=j_C2H4/self.F/12*3600*1e3/self.Mirror_factor #Convert towards umol/h/cm^2
        
        #Plot results
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        if self.plot_curves:
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,2,4,6,8,10,12,14,16],cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_Faradaic_efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$FE_{\rm C_2H_4}$: ',title='',unit='%',levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75],cmap='Blues_r',clabel=r'$FE_{\rm C_2H_4}$ (%)')
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_V_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$U_{\rm int}$: ',title='',unit='V',levels=np.arange(np.min(self.PV_EC_V_matrix),np.max(self.PV_EC_V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',clabel=r'$U_{\rm int}$ (V)',decimals=2)
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_J_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$j_{\rm s}$: ',title='',unit=' mA/cm$^2$',cmap='autumn',levels=np.array([0,3,6,9,12,15,18,21]),clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_umol_C2H4_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
            # self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_scaling_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$A_{\rm s}/A_{\rm CO_2RR}$: ',unit='x', title='',cmap='Blues_r',levels=[0,5,10,15,20,25,30,35],clabel=r'PV-EC scaling $A_{\rm s}/A_{\rm CO_2RR}$',decimals=1)
            
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Efficiency_matrix[:,:,0]*electrolyzer_efficiency*Faradaic_efficiency,xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],detail='%1.1f',cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.umol_cm2[:,:,0],xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,250/self.Solar_J_matrix[:,:,0],xmax=self.Max_efficiency_Solar_cell_2_index,ymax=self.Max_efficiency_Solar_cell_1_index,boxtext=r'$A_{\rm s}/A_{\rm CO_2RR}$: ',unit='x', title='',cmap='Blues_r',levels=[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35],clabel=r'PV-EC scaling $A_{\rm s}/A_{\rm CO_2RR}$',decimals=1)
                 
    
    def Calculate_1on1(self):
        #Calculations when no concentration is used at all
        
        self.Load_Spectrum_AM15()
        self.PEC_mode='1:1'
        self.one_one_Efficiency_matrix,self.one_one_Faradaic_efficiency_matrix,self.one_one_V_matrix,self.one_one_J_matrix,self.one_one_umol_C2H4_matrix,extra=self.Calculate_solar_cell_grid()
        if self.plot_curves:
            self.Plot_1_1_matrix()
        #
        
    def Calculate_sparse_coverage(self):
        #We assume that the Co2RR catalyst is applied perfectly on the solar cells to achieve optimal faradaic efficiency
        self.Load_Spectrum_AM15()
        self.PEC_mode='sparse coverage'
        
        self.Sparse_Coverage_Efficiency_matrix,self.Sparse_Coverage_Faradaic_efficiency_matrix,self.Sparse_Coverage_V_matrix,self.Sparse_Coverage_J_matrix,self.Sparse_Coverage_umol_C2H4_matrix,self.Sparse_Coverage_CO2_catalyst_concentration_matrix=self.Calculate_solar_cell_grid()
        if self.plot_curves:
            self.Plot_sparse_coverage_matrix()
        
        
    def Calculate_solar_light_concentration(self):
        #We assume exactly enough mirrors are used to achieve optimal faradaic efficiency
        self.Load_Spectrum_AM15()
        self.PEC_mode='solar concentration'
        self.solar_light_concentration_Efficiency_matrix,self.solar_light_concentration_Faradaic_efficiency_matrix,self.solar_light_concentration_V_matrix,self.solar_light_concentration_J_matrix,self.solar_light_concentration_umol_C2H4_matrix,self.solar_light_concentration_mirror_matrix=self.Calculate_solar_cell_grid()
        if self.plot_curves:
            self.Plot_solar_concentration_matrix()
        
    def Calculate_PV_EC(self):
        #We assume that area's between PV and EC are scaled perfectly. Here, one PV panel is used for 1 EC module. 
        
        self.Load_Spectrum_AM15()
        self.PEC_mode='PV-EC'
        self.PV_EC_Efficiency_matrix,self.PV_EC_Faradaic_efficiency_matrix,self.PV_EC_V_matrix,self.PV_EC_J_matrix,self.PV_EC_umol_C2H4_matrix,self.PV_EC_scaling_matrix=self.Calculate_solar_cell_grid()
        if self.plot_curves:
            self.Plot_PV_EC_matrix()
       
        
    
    def Calculate_solar_cell_grid(self):
        #Here, the calculation of each individual solar cell is calculated and put in a matrix, using the earlier applied settings
        
        time_start=time.time()
        
        #Set up matrices
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        if self.Three_solar_cells:
            Solar_cell_3_range=np.arange(self.bandgap_range[4],self.bandgap_range[5],self.grid_step[2])
        else:
            Solar_cell_3_range=[0]
        
        Efficiency_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        Faradaic_efficiency_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        V_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        J_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        umol_C2H4_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        if self.PEC_mode=='sparse coverage':
            CO2_catalyst_concentration_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        if self.PEC_mode=='solar concentration':
            Mirror_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
        if self.PEC_mode=='PV-EC':
            Scaling_matrix=np.zeros([np.size(Solar_cell_1_range),np.size(Solar_cell_2_range),np.size(Solar_cell_3_range)])
    
        if self.Bandgap_message==True:
            if self.Three_solar_cells:
                print('Bandgap1 \t Bandgap2 \t Bandgap3 \t STE \t umol C2H4/h/cm2 \t FE_all \t V \t\t J')
            else:
                print('Bandgap1 \t Bandgap2  \t STE \t umol C2H4/h/cm2 \t FE_all \t V \t\t J')

        
        
        
        #Calculate overlap point for every solar cell
        count_x=0
        for Bandgap_1 in Solar_cell_1_range:
            count_y=0
            for Bandgap_2 in Solar_cell_2_range:
                count_z=0
                for Bandgap_3 in Solar_cell_3_range:
                
                    #The solar cell with the bigger bandgap should go on top
                    if Bandgap_2>=Bandgap_1 or Bandgap_3>Bandgap_1 or Bandgap_3>Bandgap_2:
                        continue
                    
                    if self.PEC_mode=='solar concentration':
                        #Solar concentration is optimized to achieve the desired overpotential on the CO2RR catalyst. Not possible for a triple junction
                        self.Mirror_factor=self.J_max_ethylene/50 #We assume that no PEC cell will have a generated current density over 50 mA/cm2
                        self.Solar_Concentration=self.Mirror_factor
                        
                        V_op,J_op,Faradaic_efficiency=self.Find_current_match(Bandgap_1,Bandgap_2)
                        
                        while True:
                            if self.message_concentration:
                                print(J_op,V_op,Faradaic_efficiency,self.Mirror_factor)
                            if J_op<self.J_max_ethylene*0.9:
                                #Speed up calculations
                                self.Mirror_factor=np.ceil(self.J_max_ethylene/(J_op/self.Mirror_factor))
                            else:    
                                #Slow increase when necessary
                                self.Mirror_factor+=1
                            self.Solar_Concentration=self.Mirror_factor
                            self.Load_Spectrum_AM15()
                            V_op,J_op,Faradaic_efficiency=self.Find_current_match(Bandgap_1,Bandgap_2)
                            
                            #If less than 1 mA/cm2 is used, there is no need to calculate the ideal solar cell combinations anymore
                            if  self.Mirror_factor>self.J_max_ethylene/1:
                                self.Mirror_factor=self.J_max_ethylene/1
                                
                                V_op,J_op,Faradaic_efficiency=self.Find_current_match(Bandgap_1,Bandgap_2)
                                break
                            #This indicates that we have sufficient current density
                            if J_op>self.J_max_ethylene:
                                J_op=self.J_max_ethylene
                                Faradaic_efficiency=self.Max_ethylene_percentage
                                break
                        if self.message_concentration:
                            print(J_op,V_op,Faradaic_efficiency,self.Mirror_factor)                
                            
                            
                            
    
                    else:
                        #Find cross points between the solar cell graph and the catalyst graphs
                        if self.Three_solar_cells:
                            V_op,J_op,Faradaic_efficiency=self.Find_current_match_triple_cell(Bandgap_1,Bandgap_2,Bandgap_3)
                        
                        else:
                                
                            V_op,J_op,Faradaic_efficiency=self.Find_current_match(Bandgap_1,Bandgap_2)
                    
                    
                    #Store data
                    if self.Reduction_mode=='CO2 reduction':
                       
                        Efficiency_matrix[count_x,count_y,count_z]=self.STE(J_op,Faradaic_efficiency)
                        
                        umol_C2H4_matrix[count_x,count_y,count_z]=self.mol_ethylene(J_op,Faradaic_efficiency)
                        #print(self.Co2RcatalystConcentrator)
                        if self.PEC_mode=='sparse coverage':
                            if J_op>0:
                                CO2_catalyst_concentration_matrix[count_x,count_y,count_z]=self.J_max_ethylene/J_op
                        elif self.PEC_mode=='solar concentration':
                            Mirror_matrix[count_x,count_y,count_z]=self.Mirror_factor
                        elif self.PEC_mode=='PV-EC':
                            if J_op==0:
                                J_op=0.01
                            Scaling_matrix[count_x,count_y,count_z]=self.J_max_ethylene/J_op
                            
                            
                    elif self.Reduction_mode=='HER':
                        Efficiency_matrix[count_x,count_y,count_z]=self.STH(J_op)
                    
                    Faradaic_efficiency_matrix[count_x,count_y,count_z]=Faradaic_efficiency
                    V_matrix[count_x,count_y,count_z]=V_op
                    J_matrix[count_x,count_y,count_z]=J_op
                    
                    
                    if self.Bandgap_message==True:
                        if self.Three_solar_cells:
                
                            print('{0:.2f}'.format(Bandgap_1)+'\t\t{0:.2f}'.format(Bandgap_2)+'\t\t{0:.2f}'.format(Bandgap_3)+'\t\t\t{0:.1f}'.format(Efficiency_matrix[count_x,count_y,count_z])+'\t\t\t{0:.2f}'.format(umol_C2H4_matrix[count_x,count_y,count_z])+'\t \t{0:.2f}'.format(Faradaic_efficiency_matrix[count_x,count_y,count_z])+'\t\t{0:.2f}'.format(V_matrix[count_x,count_y,count_z])+'\t{0:.2f}'.format(J_matrix[count_x,count_y,count_z]))
                        else:
                            print('{0:.2f}'.format(Bandgap_1)+'\t\t{0:.2f}'.format(Bandgap_2)+'\t\t{0:.1f}'.format(Efficiency_matrix[count_x,count_y,count_z])+'\t\t\t\t{0:.2f}'.format(umol_C2H4_matrix[count_x,count_y,count_z])+'\t\t\t{0:.2f}'.format(Faradaic_efficiency_matrix[count_x,count_y,count_z])+'\t{0:.2f}'.format(V_matrix[count_x,count_y,count_z])+'\t{0:.2f}'.format(J_matrix[count_x,count_y,count_z]))
                      
                    count_z+=1
                count_y+=1
                
            if self.Eg_completed:    
                print('Time so far:' +'{0:.2f}'.format((time.time()-time_start)/60)+'Minutes')
            count_x+=1
        #Store data    
        self.Max_efficiency_Solar_cell_3_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[2][0]    
        self.Max_efficiency_Solar_cell_2_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[1][0]
        self.Max_efficiency_Solar_cell_1_index=np.where(Efficiency_matrix==np.max(Efficiency_matrix))[0][0]
        self.Solar_cell_combination=[Solar_cell_1_range[self.Max_efficiency_Solar_cell_1_index],Solar_cell_2_range[self.Max_efficiency_Solar_cell_2_index],Solar_cell_3_range[self.Max_efficiency_Solar_cell_3_index]]
        
        if not self.Three_solar_cells:
            #Easier for plotting
            Efficiency_matrix=Efficiency_matrix[:,:,0]
            Faradaic_efficiency_matrix=Faradaic_efficiency_matrix[:,:,0]
            V_matrix=V_matrix[:,:,0]
            J_matrix=J_matrix[:,:,0]
            umol_C2H4_matrix=umol_C2H4_matrix[:,:,0]
            if self.PEC_mode=='sparse coverage':
                CO2_catalyst_concentration_matrix=CO2_catalyst_concentration_matrix[:,:,0]
            elif self.PEC_mode=='solar concentration':
                Mirror_matrix=Mirror_matrix[:,:,0]
            elif self.PEC_mode=='PV-EC':
                Scaling_matrix=Scaling_matrix[:,:,0]
                
        
        if self.PEC_mode=='sparse coverage':
            if self.Three_solar_cells:
                self.Co2RcatalystConcentrator=CO2_catalyst_concentration_matrix[self.Max_efficiency_Solar_cell_1_index,self.Max_efficiency_Solar_cell_2_index,self.Max_efficiency_Solar_cell_3_index]
            else:
                self.Co2RcatalystConcentrator=CO2_catalyst_concentration_matrix[self.Max_efficiency_Solar_cell_1_index,self.Max_efficiency_Solar_cell_2_index]
          
            return Efficiency_matrix,Faradaic_efficiency_matrix,V_matrix,J_matrix,umol_C2H4_matrix,CO2_catalyst_concentration_matrix
        elif self.PEC_mode=='solar concentration':
            return Efficiency_matrix,Faradaic_efficiency_matrix,V_matrix,J_matrix,umol_C2H4_matrix,Mirror_matrix
        elif self.PEC_mode=='PV-EC':
            return Efficiency_matrix,Faradaic_efficiency_matrix,V_matrix,J_matrix,umol_C2H4_matrix,Scaling_matrix
        else:
            return Efficiency_matrix,Faradaic_efficiency_matrix,V_matrix,J_matrix,umol_C2H4_matrix,0
        
    def Load_Spectrum_AM15(self):
        #The photons as found in the solar spectrum, per wavelength. 
        #Source: https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html
        
        data=np.loadtxt(self.Solar_spectrum_location ,delimiter=';',dtype=str)
        data=np.char.replace(data, ',', '.')
        data = np.char.replace(data, '\'', '')
        data= np.char.replace(data,'"','')
        data = np.char.replace(data, 'b', '').astype(np.float64)
        
        #AM1.5G is in column 2, AM1.5D in column 3
        if self.Mirror_factor==1:
            self.Irradiance=data[:,2] # W/m^2/nm
        else:
            self.Irradiance=data[:,3]
        
        self.Wavelength=data[:,0]# in nm
            
        self.Irradiance=self.Solar_Concentration*self.Irradiance    #Multiply by concentration, required for variations in solar light 
        self.Spectrum_Loaded=True
        
        self.Solar_Power         = np.trapz(self.Irradiance, x=self.Wavelength)/10  # for AM1.5 solar spectrum, the total power is close to 100 mW/cm2
        self.Energy             = self.nmeV / self.Wavelength                   # eV
        self.Flux      = self.Irradiance * self.Wavelength * 1e-9 / self.h/self.c
        
        
        self.WavelengthMin      = self.Wavelength[0]
        self.WavelengthMax      = self.Wavelength[- 1]
        self.BandgapMin         = self.nmeV / self.WavelengthMax    # in eV
        self.BandgapMax         = self.nmeV / self.WavelengthMin    # in eV
        
        
        self.Second_solar_cell=False
        if self.message_spectrum_loaded:
            print('Spectrum AM1.5 loaded')
        if self.plot_spectrum:
                self.fig_spectrum=plt.figure()
                ax=self.fig_spectrum.add_subplot(1, 1, 1)
                ax.plot(self.Wavelength,self.Irradiance)
                ax.set_xlabel('Wavelength (nm)')
                ax.set_ylabel('Irradiance (W/m$^2$/nm)')
                ax.grid()
                self.fig_spectrum.tight_layout()
                #self.fig_spectrum.show()
                
                self.fig_eVvsFlux=plt.figure()
                ax=self.fig_eVvsFlux.add_subplot(1, 1, 1)
                ax.plot(self.Energy,self.Irradiance/self.Energy)
                ax.set_xlabel('Energy (eV)')
                ax.set_ylabel('Photon flux (# of photons/m$^2$/eV)')
                ax.grid()
                self.fig_eVvsFlux.tight_layout()
                #self.fig_eVvsFlux.show()
        
    
    def Find_current_match(self,Bandgap_1,Bandgap_2):
        #For a single combination of solar cells, the cross point with the earlier set catalyst is chosen
        
        #Combine solar cell data
        Voltage_SC1,J_SC1,FF,Voltage_SC2,J_SC2,FF=self.Double_solar_cell(Bandgap_1,Bandgap_2)
        V_tandem_solar_cell,J_tandem_solar_cell=self.find_J_V_tandem_solar_cel(Voltage_SC1,J_SC1,Voltage_SC2,J_SC2)
        
        
        J_catalyst=np.arange(0.1,max(J_tandem_solar_cell)*1.2,self.Current_accuracy*self.Mirror_factor)
        
        #set up catalyst graph
        if self.PEC_mode=='PV-EC':
            V_required,V_OER,V_CO2RR,Ethylene_percentage=self.Find_V_required_PEC_cell(self.J_max_ethylene)
            V_tot=np.zeros(len(J_catalyst))+V_required
            
        else:
           
            V_tot=np.zeros(len(J_catalyst))
            
            V_OER_catalyst=np.zeros(len(J_catalyst))
            V_CO2RR_catalyst=np.zeros(len(J_catalyst))
            Ethylene_percentage=np.zeros(len(J_catalyst))
            
            for i in range(len(J_catalyst)):
                J=J_catalyst[i]
                V_tot[i],V_OER_catalyst[i],V_CO2RR_catalyst[i],Ethylene_percentage[i]=self.Find_V_required_PEC_cell(J)
        #Find the intersection    
        V_op,J_op=intersect.intersection(V_tandem_solar_cell,J_tandem_solar_cell,V_tot,J_catalyst)
        
        if self.PEC_mode=='PV-EC' or self.Fixed_ethylene_output:
            Faradaic_efficiency=self.Max_ethylene_percentage
        else:
            #If sparse coverage is chose, the ethylene percentage is also assumed to be maximum at AM1.5G. Only at variations in sunlight, a interpolation is chosen
            #Faradaic_efficiency=np.interp(J_op,J_catalyst,Ethylene_percentage)
            Faradaic_efficiency=np.interp(J_op*self.Co2RcatalystConcentrator,self.Current_density_CO2RR_extended,self.Ethylene_percentage_extended)[0]
        
        if np.size(V_op)==0:
            V_op=J_op=Faradaic_efficiency=0
        
        if self.plot_intersection==True:
            plt.figure()
            plt.plot(V_tandem_solar_cell,J_tandem_solar_cell,V_tot,J_catalyst,V_op,J_op,'.k') 
        return V_op,J_op,Faradaic_efficiency
    
    def Find_current_match_triple_cell(self,Bandgap_1,Bandgap_2,Bandgap_3):
        #Similar to the code above, but with three solar cells instead of 2. 
        Voltage_SC1,J_SC1,FF,Voltage_SC2,J_SC2,FF,Voltage_SC3,J_SC3,FF=self.Triple_solar_cell(Bandgap_1,Bandgap_2,Bandgap_3)
        V_tandem_solar_cell,J_tandem_solar_cell=self.find_J_V_triple_solar_cel(Voltage_SC1,J_SC1,Voltage_SC2,J_SC2,Voltage_SC3,J_SC3)
        
        
        J_catalyst=np.arange(0.1,max(J_tandem_solar_cell)*1.2,self.Current_accuracy*self.Mirror_factor)
        if self.PEC_mode=='PV-EC':
            V_required,V_OER,V_CO2RR,Ethylene_percentage=self.Find_V_required_PEC_cell(self.J_max_ethylene)
   
            V_tot=np.zeros(len(J_catalyst))+V_required
            
        else:
           
            V_tot=np.zeros(len(J_catalyst))
            
            V_OER_catalyst=np.zeros(len(J_catalyst))
            V_CO2RR_catalyst=np.zeros(len(J_catalyst))
            Ethylene_percentage=np.zeros(len(J_catalyst))
            
            for i in range(len(J_catalyst)):
                J=J_catalyst[i]
                V_tot[i],V_OER_catalyst[i],V_CO2RR_catalyst[i],Ethylene_percentage[i]=self.Find_V_required_PEC_cell(J)
            
        V_op,J_op=intersect.intersection(V_tandem_solar_cell,J_tandem_solar_cell,V_tot,J_catalyst)

        if self.PEC_mode=='PV-EC'  or self.Fixed_ethylene_output:
            Faradaic_efficiency=self.Max_ethylene_percentage
        else:
            #Faradaic_efficiency=np.interp(J_op,J_catalyst,Ethylene_percentage)
            Faradaic_efficiency=np.interp(J_op*self.Co2RcatalystConcentrator,self.Current_density_CO2RR_extended,self.Ethylene_percentage_extended)[0]
        
        if np.size(V_op)==0:
            V_op=J_op=Faradaic_efficiency=0
        
        if self.plot_intersection:
            plt.figure()
            plt.plot(V_tandem_solar_cell,J_tandem_solar_cell,V_tot,J_catalyst,V_op,J_op,'.k') 
        return V_op,J_op,Faradaic_efficiency
        
    def Double_solar_cell(self,Bandgap_1,Bandgap_2):
        #Ensures that the second solar cell cannot abosrb photons absorbed by the first solar cell
        self.Load_Spectrum_AM15()
        self.Second_solar_cell=False
        Voltage_SC1,J_SC1,Rs1,FF1=self.Match_FF_goal(Bandgap_1)
        
        
        if Bandgap_2>=Bandgap_1: #Bigger bandgap should be on top
            return Voltage_SC1,J_SC1,Voltage_SC1,0*Voltage_SC1
        
        self.Second_solar_cell=True
        Voltage_SC2,J_SC2,Rs2,FF2=self.Match_FF_goal(Bandgap_2)
        
        return Voltage_SC1,J_SC1,FF1,Voltage_SC2,J_SC2,FF2
    
    def Triple_solar_cell(self,Bandgap_1,Bandgap_2,Bandgap_3):
        #Ensures that the second (and third) solar cell cannot abosrb photons absorbed by the first (and second) solar cell
        self.Load_Spectrum_AM15()
        self.Second_solar_cell=False
        Voltage_SC1,J_SC1,Rs1,FF1=self.Match_FF_goal(Bandgap_1)
        
        
        if Bandgap_2>=Bandgap_1 or Bandgap_3>=Bandgap_1: #Bigger bandgap should be on top
            return Voltage_SC1,J_SC1,FF1,Voltage_SC1,0*Voltage_SC1,0*FF1,Voltage_SC1,0*Voltage_SC1,0*Voltage_SC1
        
        self.Second_solar_cell=True
        Voltage_SC2,J_SC2,Rs2,FF2=self.Match_FF_goal(Bandgap_2)
        
        if Bandgap_3>=Bandgap_2: #Bigger bandgap should be on top
            return Voltage_SC1,J_SC1,FF1,Voltage_SC2,J_SC2,FF2,Voltage_SC1,0*Voltage_SC1,0*Voltage_SC1
        
        
        #Set spectrum to 0 after the second solar cell
        self.Second_solar_cell=False
        Rubish=self.Solar_cell(Bandgap_2)
        
        self.Second_solar_cell=True
        Voltage_SC3,J_SC3,Rs3,FF3=self.Match_FF_goal(Bandgap_3)
        return Voltage_SC1,J_SC1,FF1,Voltage_SC2,J_SC2,FF2,Voltage_SC3,J_SC3,FF3
    
    def Match_FF_goal(self,Bandgap):
        #Using Rs, the solar cell graph is adjusted untill the desired FF is found. 
        
        V,J=self.Solar_cell(Bandgap,Rs=0)
        FF=self.FF_Solar_cell(V,J)
        Voc_ideal=intersect.intersection(V,J,V, V*0)[0][0]
        
        if FF<self.FF_goal:
            return V,J,0,FF
        #We start slow, but the Rs increases rapidly. 
        Rs_start=1e-10
        Rs=Rs_start
        if Rs<1e-10:
            Rs=1e-10
        previous='greater'
        count_FF=0
        R_factor=4
        #This part could probably be more optimized if the code needs to speed up. 
        while abs(FF-self.FF_goal)>0.001 and count_FF<100:
            
            if FF>self.FF_goal:
                if previous=='smaller':
                    R_factor=1+(R_factor-1)/2
                Rs*=R_factor
                previous='greater'
            elif FF<self.FF_goal:
                if previous=='greater':
                    R_factor=1+(R_factor-1)/2
                Rs/=R_factor
                previous='smaller'
            
            J=self.Solar_cell_total_current_Lambert_W(Bandgap,V,J[0],Rs)
            FF=self.FF_Solar_cell(V, J)
         
            if Rs<1e-12:
                Rs=1e-12
                break
            count_FF+=1
        Voc_end=intersect.intersection(V,J,V, V*0)[0][0]
        if abs(Voc_end-Voc_ideal)>1e-2:
            print('Warning, Voc too low')
       
        return V,J,Rs,FF
        
        
    def Solar_cell(self,Bandgap,Rs=0,Rsh=1e99,Bandgap_Cutoff = 0.0):
        #Calculates a individual solar cell graph
        
        Lambda_Max         = self.nmeV / Bandgap
        Lambda_Min     = 0.0
        if Bandgap_Cutoff>0:
            Lambda_Min  = self.nmeV / Bandgap_Cutoff
        # end if
        
        indexT         = np.where(np.logical_and(self.Wavelength >= Lambda_Min, self.Wavelength <= Lambda_Max))
        Wavelength     = self.Wavelength[indexT]
        #Energy         = self.Energy[indexT]
        
            
        if self.Second_solar_cell:
            #Here we already have a second or third solar cell
            Flux=self.Flux_after[indexT]
            
            #These have different reflections as well
            self.n_top=3.3 #Value for GaAs
            self.n_bottom=0 #perfect back reflector
        else:
            #Set the spectrum for the second (or third) solar cell  
            Flux           = self.Flux[indexT]
            self.Flux_after=copy.copy(self.Flux)
            self.Flux_after[indexT]=0
            self.n_top=1 #Value for Si
            self.n_bottom=3.42 #air
        
        Jph           = self.q * np.trapz(Flux, x=Wavelength) /10      # Short-Circuit Current in mA/cm2
        J0_ideal      =self.J_rad_Henri_diode(Bandgap,0,0)

        
        VOC_ideal            = self.kTeV * np.log((Jph / J0_ideal) + 1.0)         # Open-Circuit Voltage in V
        Vstep          = max(VOC_ideal / 500.0,0.0001)
        Voltage        = np.arange(0.0, VOC_ideal + 2*Vstep, Vstep)

        

            
        J=self.Solar_cell_total_current_Lambert_W(Bandgap,Voltage,Jph,Rs,Rsh)

        return (Voltage, J)
    
    def J_rad_Henri_diode(self,Bandgap,Voltage,J,Rs=0,A_diode=1):
        #Calculate J_rad, used to determinate j0
        #C. H. Henry, Journal of Applied Physics 51, 4494 (1980); https://doi.org/10.1063/1.328272
        A=(self.q*(self.n_top**2+self.n_bottom**2)*(Bandgap*self.q)**2*self.kb*self.Temperature)/(4*np.pi**2*(self.h/2/np.pi)**3*self.c**2)/10 #A/cm^2
        exponent=np.exp((self.q*(Voltage+J*Rs-Bandgap))/(A_diode*self.kb*self.Temperature))
        return A*exponent
    
    def Solar_cell_total_current_Lambert_W(self,Bandgap,V,Jph,Rs=0,Rsh=1e99,A_diode=1,W=0.99,label=''):
        #Calculate the total response current, based on Jph and J0
        
        Vt=self.kb*self.Temperature/self.q
        J0=self.J_rad_Henri_diode(Bandgap,0,0)
     
        if Rs<=0:
            return Jph-J0*(np.exp(V/A_diode/Vt)-1)-V/Rsh
        else:
            argW=Rs*J0/(A_diode*Vt*(1+Rs/Rsh))*np.exp(V/A_diode/Vt*(1-Rs/(Rs+Rsh))+((Jph+J0)*Rs)/(A_diode*Vt*(1+Rs/Rsh)))
            return (Jph+J0-V/Rsh)/(1+Rs/Rsh)-A_diode*Vt/Rs*lambertw(argW).real
        
    def FF_Solar_cell(self,Voltage,Current):
        #Calculate the fill factor
        Pm     = Current*Voltage
        try:
            VOC=intersect.intersection(Voltage, Current, Voltage, Voltage*0)[0][0]
        except:
            VOC=np.max(Voltage)
        #Vm     = Voltage[np.where(Pm==np.max(Pm))]
        #Jm     = Current[np.where(Pm==np.max(Pm))]
        FF     = np.max(Pm) / (np.max(Current) * VOC)
        return FF
    
    def Efficiency_solar_cell(self,Voltage,Current):
        #Maximum solar to electricity efficiency of a solar cell
        Pm=Current*Voltage
        index=np.where(Pm==max(Pm))[0][0]

        return np.max(Pm)/self.Solar_Power,Voltage[index],Current[index]
        
        
    
    def find_J_V_tandem_solar_cel(self,V1,J1,V2,J2):
        #Combining two solar cell
        length=length=min(len(J1),len(J2))
        final_J=np.zeros(length)
        final_V=np.zeros(length)
        
        for i in range(length): 
            final_J[i]=J1[i] #max_current SC1 is always lower than SC2
            final_V[i]=V1[i]+np.interp(J1[i],J2[::-1],V2[::-1])
        
        return final_V,final_J
        
    def find_J_V_triple_solar_cel(self,V1,J1,V2,J2,V3,J3):
        #Similar to two solar cells, but now with three
        V23,J23=self.find_J_V_tandem_solar_cel(V2, J2, V3, J3)
        
        final_V,final_J=self.find_J_V_tandem_solar_cel(V1, J1, V23, J23)     

        return final_V,final_J
    
    
        
        
    def Find_V_required_PEC_cell(self,J):
        #Data for the V,J plot for catalysts, based on impunt current. 
        
        V_OER_catalyst=self.Voltage_OER(J)
        
        if self.Reduction_mode=='CO2 reduction':
            #Normally, we would assume that the design is so that the faradaic efficiency is maximal for the CO2RR catalyst. However, for changes in light conditions, the faradaic efficiency can vary. 
            if self.PEC_mode=='sparse coverage' and self.Fixed_ethylene_output:
                
                V_CO2RR_catalyst,Ethylene_percentage=self.V_CO2RR(J,optimal=True)
            else:
                
                V_CO2RR_catalyst,Ethylene_percentage=self.V_CO2RR(J*self.Co2RcatalystConcentrator,optimal=False)
            V_tot=V_OER_catalyst+V_CO2RR_catalyst+J*1e-3*self.Fluid_resistance
            return V_tot,V_OER_catalyst,V_CO2RR_catalyst,Ethylene_percentage
        
        
        elif self.Reduction_mode=='HER':
            V_HER_catalyst=self.V_HER(J)
            V_tot=V_OER_catalyst+V_HER_catalyst+J*1e-3*self.Fluid_resistance
            return V_tot,V_OER_catalyst,V_HER_catalyst,100
            
    def Voltage_OER(self,J):
        
        eta_OER=self.Tafel_slope(self.i_OER,self.eta_OER,self.A_OER,J)
        
        return 1.23+eta_OER
    
    
    
    def V_CO2RR(self,J,optimal=False):
        
        if optimal:
            #Fixed overpotential & ethylene percentage
            return -self.V_max_ethylene,self.Max_ethylene_percentage
        
        else:
            #Find the linear interpolation
            
            location=np.where(abs(self.Current_density_CO2RR_extended-J)==min(abs(self.Current_density_CO2RR_extended-J)))[0][0]
            
            return -self.eta_CO2RR[location],self.Ethylene_percentage_extended[location]
            
            
            
    def V_HER(self,J):
        return self.Tafel_slope(self.i_HER,self.eta_HER,self.A_HER,J)
        
        
    def Tafel_slope(self,i_1,eta_1,A,J)  :
        #Use only when eta>0.1. Give a standard value for i0 and eta, the tafel slope and your current and it will return the needed overpotential in V)
        i0=i_1/(10**(eta_1/A))
        return A*np.log10(J/i0)
    
        
        
        
        
        
        
        
        
    
    def STE(self,J_op,Faradaic_efficiency):
        #Solar to Ethylene Efficiency
        return J_op*(self.V_OER-self.V_CO2R)/self.Solar_Power*Faradaic_efficiency
    
    def mol_ethylene(self,J_op,Faradaic_efficiency):
        #Solar to ethylene production in umol/h/cm2
        J_C2H4=J_op*Faradaic_efficiency/100 #mA/cm^2 towards etylene
        return J_C2H4/self.F/12*3600*1e3/self.Mirror_factor #Convert towards umol/h/cm^2
    
    def STH(self,J_op):
        #Solar to hydrogen efficiency
        return J_op*(self.V_OER-self.V_HER_potential)/self.Solar_Power*100
        
    
    def Plot_1_1_matrix(self):
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        Max_efficiency_Solar_cell_2=np.where(self.one_one_Efficiency_matrix==np.max(self.one_one_Efficiency_matrix))[1][0]
        Max_efficiency_Solar_cell_1=np.where(self.one_one_Efficiency_matrix==np.max(self.one_one_Efficiency_matrix))[0][0]
    
    
        if self.Reduction_mode=='HER': 
          self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext='no',unit='%',levels=[0.5,3.5,6.5,9.5,12.5,15.5,18.5,21.5,24.5,27.5,30.5],cmap='Greens_r',clabel=r'$\eta_{\rm STH}$ (%)')
        else:
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,2,4,6,8,10,12,14,16],cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_Faradaic_efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$FE_{\rm C_2H_4}$: ',title='',unit='%',levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75],cmap='Blues_r',clabel=r'$FE_{\rm C_2H_4}$ (%)')
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_V_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$U_{\rm int}$: ',title='',unit='V',levels=np.arange(np.min(self.one_one_V_matrix),np.max(self.one_one_V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',clabel=r'$U_{\rm int}$ (V)',decimals=2)
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_J_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$j_{\rm s}$: ',title='',unit=' mA/cm$^2$',cmap='autumn',levels=[0,3,6,9,12,15,18,21],clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
            self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.one_one_umol_C2H4_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
            
        
    
    def Plot_sparse_coverage_matrix(self):
        
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        Max_efficiency_Solar_cell_2=np.where(self.Sparse_Coverage_Efficiency_matrix==np.max(self.Sparse_Coverage_Efficiency_matrix))[1][0]
        Max_efficiency_Solar_cell_1=np.where(self.Sparse_Coverage_Efficiency_matrix==np.max(self.Sparse_Coverage_Efficiency_matrix))[0][0]
        
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,2,4,6,8,10,12,14,16],cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_Faradaic_efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$FE_{\rm C_2H_4}$: ',title='',unit='%',levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75],cmap='Blues_r',clabel=r'$FE_{\rm C_2H_4}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_V_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$U_{\rm int}$: ',title='',unit='V',levels=np.arange(np.min(self.Sparse_Coverage_V_matrix),np.max(self.Sparse_Coverage_V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',decimals=2,clabel=r'$U_{\rm int}$ (V)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_J_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$j_{\rm s}$: ',title='',unit=' mA/cm$^2$',cmap='autumn',levels=[0,3,6,9,12,15,18,21],clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_umol_C2H4_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.Sparse_Coverage_CO2_catalyst_concentration_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$A_{\rm s}/A_{\rm CO_2RR}$: ',unit='', title='',cmap='Reds',levels=[0,5,10,15,20,25,30,35,40,45,50],clabel=r'$A_{\rm s}/A_{\rm CO_2RR}$',decimals=1)
        
        
    def Plot_solar_concentration_matrix(self):
        
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        Max_efficiency_Solar_cell_2=np.where(self.solar_light_concentration_Efficiency_matrix==np.max(self.solar_light_concentration_Efficiency_matrix))[1][0]
        Max_efficiency_Solar_cell_1=np.where(self.solar_light_concentration_Efficiency_matrix==np.max(self.solar_light_concentration_Efficiency_matrix))[0][0]
        
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,2,4,6,8,10,12,14,16],cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_Faradaic_efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$FE_{\rm C_2H_4}$: ',title='',unit='%',levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75],cmap='Blues_r',clabel=r'$FE_{\rm C_2H_4}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_V_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$U_{\rm int}$: ',title='',unit='V',levels=np.arange(np.min(self.solar_light_concentration_V_matrix),np.max(self.solar_light_concentration_V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',clabel=r'$U_{\rm int}$ (V)',decimals=2)
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_J_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$j_{\rm s}$: ',title='',unit=' mA/cm$^2$',cmap='autumn',levels=np.array([0,3,6,9,12,15,18,21])*self.J_max_ethylene/20,clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_umol_C2H4_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.solar_light_concentration_mirror_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext='$C$: ',unit='', title='',cmap='Blues_r',levels=[0,5,10,15,20,25,30,35],clabel='Solar light concentration $C$',decimals=1)
           
    def Plot_PV_EC_matrix(self):
    
        Solar_cell_1_range=np.arange(self.bandgap_range[0],self.bandgap_range[1],self.grid_step[0])
        Solar_cell_2_range=np.arange(self.bandgap_range[2],self.bandgap_range[3],self.grid_step[1])
        Max_efficiency_Solar_cell_2=np.where(self.PV_EC_Efficiency_matrix==np.max(self.PV_EC_Efficiency_matrix))[1][0]
        Max_efficiency_Solar_cell_1=np.where(self.PV_EC_Efficiency_matrix==np.max(self.PV_EC_Efficiency_matrix))[0][0]
        
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_Efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\eta_{\rm STE}$: ',unit='%',levels=[0,2,4,6,8,10,12,14,16],cmap='Greens_r',clabel=r'$\eta_{\rm STE}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_Faradaic_efficiency_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$FE_{\rm C_2H_4}$: ',title='',unit='%',levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75],cmap='Blues_r',clabel=r'$FE_{\rm C_2H_4}$ (%)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_V_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$U_{\rm int}$: ',title='',unit='V',levels=np.arange(np.min(self.PV_EC_V_matrix),np.max(self.PV_EC_V_matrix)+0.05,0.05),detail='%1.2f',cmap='winter',clabel=r'$U_{\rm int}$ (V)',decimals=2)
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_J_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$j_{\rm s}$: ',title='',unit=' mA/cm$^2$',cmap='autumn',levels=np.array([0,3,6,9,12,15,18,21]),clabel=r'$j_{\rm s}$ (mA/cm$^2$)',decimals=1)
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_umol_C2H4_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$\nu_{\rm C_2H_4}$: ',unit=r' $\mu$mol/h/cm$^2$', title='',cmap='Purples_r',levels=[0,4,8,12,16,20,24,28,32,36],clabel=r'$\nu_{\rm C_2H_4}$ ($\mu$mol/h/cm$^2$)')
        self.Plot_matrix(Solar_cell_1_range,Solar_cell_2_range,self.PV_EC_scaling_matrix,xmax=Max_efficiency_Solar_cell_2,ymax=Max_efficiency_Solar_cell_1,boxtext=r'$A_{\rm s}/A_{\rm CO_2RR}$: ',unit='x', title='',cmap='Blues_r',levels=[0,5,10,15,20,25,30,35],clabel=r'PV-EC scaling $A_{\rm s}/A_{\rm CO_2RR}$',decimals=1)
              
       
        
    def Plot_matrix(self,SC_1_set,SC_2_set,Z,xmax=0,ymax=0,boxtext='STE: ',unit='%',levels=[0],detail='%1.1f',cmap='brg',title='',clabel='',decimals=0):
        X=np.meshgrid(SC_2_set,SC_1_set)[0]
        Y=np.meshgrid(SC_2_set,SC_1_set)[1]
        plt.figure()
        figcontour, ax = plt.subplots()
        
        if max(levels)==0:
            levels= np.arange(np.min(Z),np.max(Z)*1.1,(np.max(Z)-np.min(Z))/10)

        CS1 = ax.contourf(X, Y, Z,levels=levels,cmap=cmap)  
        # CS1 = ax.contourf(X, Y, Z,levels=np.linspace(np.min(levels),np.max(levels),1000),cmap=cmap)
        CS2 = ax.contour(X, Y, Z,levels=levels,colors='w')
        cbar=figcontour.colorbar(CS1,ticks=levels)    
        cbar.ax.set_ylabel(clabel,fontsize=14)
        if xmax==0 and ymax==0:
            xmax=np.where(Z==np.max(Z))[1][0]
            ymax=np.where(Z==np.max(Z))[0][0]
            
        plt.plot(SC_2_set[xmax],SC_1_set[ymax],'ok')
        ax.set_title(title)
        ax.set_xlabel(r'$\it{E}$$_{\rm g2}$ (eV)',fontsize=14)
        ax.set_ylabel(r'$\it{E}$$_{\rm g1}$ (eV)',fontsize=14)
        props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=1)
        if boxtext!='no':
            if decimals==1:
                ax.text(1.05*min(SC_2_set),0.95*max(SC_1_set),'{0:.2f}'.format(SC_1_set[ymax])+'eV/'+'{0:.2f}'.format(SC_2_set[xmax])+'eV \n'+boxtext+'{0:.1f}'.format(Z[ymax,xmax])+unit,bbox=props,fontsize=12)
            elif decimals==2:
                ax.text(1.05*min(SC_2_set),0.95*max(SC_1_set),'{0:.2f}'.format(SC_1_set[ymax])+'eV/'+'{0:.2f}'.format(SC_2_set[xmax])+'eV \n'+boxtext+'{0:.2f}'.format(Z[ymax,xmax])+unit,bbox=props,fontsize=12)
            elif decimals==0:
                ax.text(1.05*min(SC_2_set),0.95*max(SC_1_set),'{0:.2f}'.format(SC_1_set[ymax])+'eV/'+'{0:.2f}'.format(SC_2_set[xmax])+'eV \n'+boxtext+'{0:.0f}'.format(Z[ymax,xmax])+unit,bbox=props,fontsize=12)
            else:
                #Default 1 decimal
                ax.text(1.05*min(SC_2_set),0.95*max(SC_1_set),'{0:.2f}'.format(SC_1_set[ymax])+'eV/'+'{0:.2f}'.format(SC_2_set[xmax])+'eV \n'+boxtext+'{0:.1f}'.format(Z[ymax,xmax])+unit,bbox=props,fontsize=12)
          
                
    def Vary_Solar_Intensity(self,Bandgap_1,Bandgap_2,Solar_Concentration_range=np.arange(0.1,1.5,0.01),plot_graph_solar=True):
        #Used for variations in solar efficiency
        self.Fixed_ethylene_output=False
        Efficiencies=np.zeros(len(Solar_Concentration_range))
        FEs=np.zeros(len(Solar_Concentration_range))
        Voltages=np.zeros(len(Solar_Concentration_range))
        Currents=np.zeros(len(Solar_Concentration_range))
        umol_C2H4=np.zeros(len(Solar_Concentration_range))
        i=0
        if self.Mirror_factor==1:
            #AM1.5G
            Standard_Solarpower=100.2
            Concentration_Range=Solar_Concentration_range
        else:
            #AM1.5D
            Standard_Solarpower=90*self.Mirror_factor
            Concentration_Range=Solar_Concentration_range*self.Mirror_factor
            
            
        for solar_concentration in Concentration_Range:

            
            self.Solar_Concentration=solar_concentration
            V_op,J_op,Faradaic_efficiency=self.Find_current_match(Bandgap_1,Bandgap_2)
            if self.message_variation:
                print(J_op*self.Co2RcatalystConcentrator,V_op,Faradaic_efficiency)
            

            Efficiencies[i]=self.STE(J_op,Faradaic_efficiency)
            FEs[i]=Faradaic_efficiency
            Voltages[i]=V_op
            Currents[i]=J_op
            
            umol_C2H4[i]=self.mol_ethylene(J_op, Faradaic_efficiency)
            
            i+=1
    
        if plot_graph_solar:
            
            self.plot_graph_vary_intensity(Solar_Concentration_range *Standard_Solarpower, Efficiencies, 'Solar insolation (mW/cm$^2$)',r'$\eta_{\rm STE}$ (%)')
            self.plot_graph_vary_intensity(Solar_Concentration_range *Standard_Solarpower, FEs, 'Solar insolation (mW/cm$^2$)',r'Overall FE_${\rm C_2H_4}$ (%)')
            self.plot_graph_vary_intensity(Solar_Concentration_range *Standard_Solarpower, Voltages, 'Solar insolation (mW/cm$^2$)',r'U$_{int}$')
            self.plot_graph_vary_intensity(Solar_Concentration_range *Standard_Solarpower, Currents, 'Solar insolation (mW/cm$^2$)','j_$s$ (mA/cm$^2$)',linear=1)
            self.plot_graph_vary_intensity(Solar_Concentration_range *Standard_Solarpower,  umol_C2H4, 'Solar insolation (mW/cm$^2$)',r'$\nu_{\rm C_2H_4}$ ($\mu$mol / h / cm$^2$)',linear=1)
            
            
    
        return Efficiencies,FEs,Voltages,Currents,umol_C2H4
    
    
    def Degradation_catalyst(self,Bandgap_1,Bandgap_2,Degradation_range=np.arange(0,36,1),plot_graph=True):
        #Used for variations in catalyst degradation, only for PEC



        self.Fixed_ethylene_output=False
        Efficiencies=np.zeros(len(Degradation_range))
        FEs=np.zeros(len(Degradation_range))
        Voltages=np.zeros(len(Degradation_range))
        Currents=np.zeros(len(Degradation_range))
        umol_C2H4=np.zeros(len(Degradation_range))
        j=0
        plt.figure(1)
        if self.Mirror_factor==1:
            #AM1.5G
            Standard_Solarpower=100.2
        else:
            #AM1.5D
            Standard_Solarpower=90*self.Mirror_factor
            
            
        for degradation_percentage in Degradation_range:
            
            Voltage_SC1,J_SC1,FF,Voltage_SC2,J_SC2,FF=self.Double_solar_cell(Bandgap_1,Bandgap_2)
            V_tandem_solar_cell,J_tandem_solar_cell=self.find_J_V_tandem_solar_cel(Voltage_SC1,J_SC1,Voltage_SC2,J_SC2)
            
            #Calculate jV OER
            #Calculate jV+FE CO2RR
            #Find match again with
        
            J_catalyst=np.arange(0.1,max(J_tandem_solar_cell)*1.2,self.Current_accuracy*self.Mirror_factor)
            V_tot=np.zeros(len(J_catalyst))
            
            V_OER_catalyst=np.zeros(len(J_catalyst))
            V_CO2RR_catalyst=np.zeros(len(J_catalyst))
            Ethylene_percentage=np.zeros(len(J_catalyst))
            new_Co2RcatalystConcentrator=self.Co2RcatalystConcentrator*(100/(100-degradation_percentage))
            for i in range(len(J_catalyst)):
                J_solar=J_catalyst[i]
                J_deg=J_solar*(100/(100-degradation_percentage)) #Transform the solar current density to the degradation current density, without taking concentration into account
                V_OER_catalyst[i]=self.Voltage_OER(J_deg)
                # V_OER_catalyst[i]=self.Voltage_OER(J_solar)
                V_CO2RR_catalyst[i],Ethylene_percentage[i]=self.V_CO2RR(J_deg*new_Co2RcatalystConcentrator,optimal=False)
                V_tot[i]=V_OER_catalyst[i]+V_CO2RR_catalyst[i]+J_solar*1e-3*self.Fluid_resistance
                
                
                
         
            #Find the intersection, with J_catalyst adjusted for degradation- more current has to pass through it. This assumes equal degradation of the OER and CO2RR catalysts
            V_op,J_op=intersect.intersection(V_tandem_solar_cell,J_tandem_solar_cell,V_tot,J_catalyst)
            if np.size(V_op)==0:
                V_op=J_op=Faradaic_efficiency=0
            else:
            
                Faradaic_efficiency=np.interp(J_op*new_Co2RcatalystConcentrator,self.Current_density_CO2RR_extended,self.Ethylene_percentage_extended)[0]
                
            if j==0:
                
                plt.plot(V_tandem_solar_cell,J_tandem_solar_cell,'k',V_tot,J_catalyst,'--k',V_op,J_op,'.k')
        
                
           
        
                
            
        
            if self.message_variation:
                print(J_op,J_op*new_Co2RcatalystConcentrator,V_op,Faradaic_efficiency)
            
        
            Efficiencies[j]=self.STE(J_op,Faradaic_efficiency)[0]
            FEs[j]=Faradaic_efficiency
            Voltages[j]=V_op[0]
            Currents[j]=J_op[0]
            umol_C2H4[j]=self.mol_ethylene(J_op, Faradaic_efficiency)
            
            j+=1
            
        
        plt.figure(1)
        plt.plot(V_tot,J_catalyst,':k',V_op,J_op,'.k')
        plt.xlabel('$U$ (V)',fontsize=14)
        plt.ylabel(r'$J/A_{\rm s}$',fontsize=14)
        plt.xlim(1.0,3.0)
        plt.show()
        
        
        if plot_graph:
            
            self.plot_graph_vary_intensity(Degradation_range, Efficiencies, 'Degradation of the catalysts (%)',r'$\eta_{\rm STE}$ (%)')
            self.plot_graph_vary_intensity(Degradation_range, FEs, 'Degradation of the catalysts (%)',r'Overall FE_${\rm C_2H_4}$ (%)')
            self.plot_graph_vary_intensity(Degradation_range, Voltages, 'Degradation of the catalysts (%)',r'U$_{\rm int}$')
            self.plot_graph_vary_intensity(Degradation_range, Currents, 'Degradation of the catalysts (%)',r'$j_{\rm s}$ (mA/cm$^2$)')
            self.plot_graph_vary_intensity(Degradation_range, Currents*100/(100-Degradation_range)*self.Co2RcatalystConcentrator, 'Degradation of the catalysts (%)',r'$j_{\rm CO_2RR}$ (mA/cm$^2$)')
            self.plot_graph_vary_intensity(Degradation_range,  umol_C2H4, 'Degradation of the catalysts (%)',r'$\nu_{\rm C_2H_4}$ ($\mu$mol / h / cm$^2$)')
            
            
            
            plt.show()
         
        return Efficiencies,FEs,Voltages,Currents,umol_C2H4
        
    def plot_graph_vary_intensity(self,x,y,xlabel='',ylabel='',title='',linear=0):
        
        
        if self.Mirror_factor==1:
            Design_value=100.2*self.Mirror_factor
        else:
            Design_value=90*self.Mirror_factor
        plt.figure()
        self.figcontour, ax = plt.subplots()
        if linear==1:
            y_value_100W=y[np.where(abs(x-Design_value)==min(abs(x-Design_value)))[0][0]]
            y_linear=np.linspace(0,y_value_100W*max(x)/Design_value,10)
            plt.plot(np.linspace(0,max(x),10),y_linear,'k',label='Linear interpolation of C$_2$H$_4$ production rate')
            
        plt.plot(x,y,'c',label='C$_2$H$_4$ production rate',)

        # ax.legend(fontsize=13)
        ax.set_xlabel(xlabel,fontsize=14)
        ax.set_ylabel(ylabel,fontsize=14)
        ax.set_title(title)  

        






 



