**Code used to calculate theoretical PEC cells, both for hydrogen and ethylene production for the paper: 
"Axiomatic Design for Analysis and Noniterative Engineering of Photoelectrochemical Cells for CO2 Conversion at High Photon to Product Yield"**

The main code is in 'PEC_cell.py'. Seven calculations are preselected. Scenario A uses sparse coverage (Figure 3), Scenario B uses solar concentration (Figure S2), Scenario C uses PV-EC (Figure S3) and scenario F no concentration (Figure S7), all calculated with tandem solar cells. Scenario T uses sparse coverage with tripple solar cells (Figure S11). All scenarios use data for the O(II)D-Cu and NiFeOx catalyst (Asiri, 2022 & McCrory, 2013). Variation in solar light and degradation of the surface area of the catalysts is presented for both of the Cu-Ag catalyst (Gao 2019) and the O(II)D-Cu catalyst (Figures 4, S5 and S6). Isolines for overpotential and Faradaic Efficiencies of the CO2RR catalyst (Figure 5 and Figure S7).  Vary FF solar cells calculates Scenario A with solar cells using FF=0.75 and FF=0.65 instead of the standard assumed FF=0.85 (Figure S9).    

Four types of concentration are distinguished: sparse coverage, solar concentration, PV-EC and no concentration, each with their own module. 
The reduction reaction is assumed to be CO2RR towards ethylene, but the model can also be used for hydrogen evolution reaction (HER) calculations. 
Different CO2RR catalysts for ethylene production are parameterized in the required voltage, current density and percentage ethylene, similar for OER catalysts.
Additionally, the fluid resistance (standard: 5Ω) and FF (standard 0.85) are specified. 

Using these parameters and the AM1.5G or AM1.5D spectrum23, V, j_s,FE_C2H4,ν_C2H4 and η_STE are calculated for different combinations of tandem solar cells. 
The high bandgap solar cell is assumed to be transparent for all wavelengths above its absorption threshold, allowing these photons to be absorbed by a medium bandgap solar cell. For each available current density, the voltages of the first and second solar cells are added together, resulting in a tandem solar performance graph. The intersection of the tandem solar cells graph with the catalyst graph gives j_s, used to calculate V, j_s,FE_C2H4 (based on interpolated data), ν_C2H4 and η_STE. For the OER catalyst, it is standard assumed j_s=j_OER, and J=j_s⋅1cm2. 

For the scenarios where the optimal solar cells are calculated with sparse coverage, the CO2RR catalyst is assumed to operate at the overpotential that generates the highest Faradaic efficiency. For every solar cell combination, a different j_s is calculated, which in turn is used to determine the concentration ratio A_s/A_(CO_2 RR) for the sparsely applied catalyst scenario. In the concentrated sunlight scenario, the solar concentration on the solar cells is gradually increased to achieve the ideal j_CO2RR=j_s, optimizing the selectivity towards ethylene. Additionally, R_transport is lowered to 1Ω to account for the larger surface area of the CO2RR catalyst. The optimal solar cells for electricity production are calculated and matched to an EC cell that operates at a fixed voltage. The scaling, both of the different area and number of PV and EC components, is determined afterwards. The efficiencies are based on the area of the solar cells. 
When varying sunlight and calculating the scenario without concentration, the Faradaic  efficiency, voltage and current responses of the CO2RR catalysts are interpolated from reported data.

**References**

**Calculations:**  
C. H. Henry, Journal of Applied Physics 51, 4494 (1980); https://doi.org/10.1063/1.328272

Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

Rühle, S. Tabulated Values of the Shockley–Queisser Limit for Single Junction Solar Cells. Solar Energy 2016, 130, 139–147; https://doi.org/10.1016/j.solener.2016.02.015.

Solar Spectrum: https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html  

Intersect:https://github.com/sukhbinder/intersection/blob/master/intersect/intersect.py

**CO2RR to ethylene catalysts:**  
Dan Ren ACS Sustainable Chemistry & Engineering,5(10):9191–9199 (2017); https://doi.org/10.1021/acssuschemeng.7b02110

Cao-Thang Dinh et al. ,CO2 electroreduction to ethylene via hydroxide-mediated copper catalysis at an abrupt interface.Science360,783-787(2018).DOI:10.1126/science.aas9100

Ren D, Gao J, Zakeeruddin SM, Grätzel M. New Insights into the Interface of Electrochemical Flow Cells for Carbon Dioxide Reduction to Ethylene. J Phys Chem Lett. 2021 Aug 12;12(31):7583-7589. doi: 10.1021/acs.jpclett.1c02043. Epub 2021 Aug 4. PMID: 34347495.

Jing Gao, Hong Zhang, Xueyi Guo, Jingshan Luo, Shaik M. Zakeeruddin, Dan Ren*, and Michael Grätzel*J. Am. Chem. Soc. 2019, 141, 47, 18704–18714

Jing Gao  Dan Ren Xueyi Guo Shaik Mohammed Zakeeruddin and  Michael Grätzel  Faraday Discuss., 2019,215, 282-296

Dongxing Tan, Bari Wulan, Xueying Cao, Jintao Zhang, Nano Energy,Volume 89, Part B, 2021, 106460, 2211-2855 https://doi.org/10.1016/j.nanoen.2021.106460.

Choi, C., Kwon, S., Cheng, T. et al. Highly active and stable stepped Cu surface for enhanced electrochemical CO2 reduction to C2H4. Nat Catal 3, 804–812 (2020). https://doi.org/10.1038/s41929-020-00504-x

Asiri, A. M.; Gao, J.; Khan, S. B.; Alamry, K. A.; Marwani, H. M.; Khan, M. S. J.; Adeosun, W. A.; Zakeeruddin, S. M.; Ren, D.; Grätzel, M. Revisiting the Impact of Morphology and Oxidation State of Cu on CO2 Reduction Using Electrochemical Flow Cell. J. Phys. Chem. Lett. 2022, 13 (1), 345–351. https://doi.org/10.1021/acs.jpclett.1c03957.

**HER catalyst:**  
Shu Hu et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

**OER catalyst:**  
Shu H et al., Energy Environ. Sci., ,6, 2984-2993 (2013); https://doi.org/10.1039/C3EE40453F

Joya & de Groot ACS Catal. 2016, 6, 3, 1768–1771

McCrory, C. C. L.; Jung, S.; Peters, J. C.; Jaramillo, T. F. Benchmarking Heterogeneous Electrocatalysts for the Oxygen Evolution Reaction. J. Am. Chem. Soc. 2013, 135 (45), 16977–16987. https://doi.org/10.1021/ja407115p.  





@author: assensme

# AD-for-a-PEC-cell
