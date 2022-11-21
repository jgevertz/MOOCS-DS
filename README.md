# MOOCS-DS
Multi-Objective Optimization of Combination Synergy - Dose Selection Algorithm
Data and code associated with "Guiding model-driven combination dose selection using multi-objective synergy optimization"
- All paper results related to the toy model are found by running MOOCS_DS_ToyModel.m in the ToyModel folder. Note you do need a version of MATLAB that can run a script for which the function files are embedded in the main m-file. If your version of MATLAB does not allow this, each function file needs to be put into its own .m file in the same directory and appropriately named.
- For the Pembro + Avastin model, runn MOOCS_DS_PembroAvastin.m in the PembroAvastin folder. This creates all plots using the experimental protocol (Q3D). Note the same comment about embedded functions applies here. - Synery_Across_Protocols.m 
- To generate and visualize results for protocols other than the experimental protocol, run Synergy_Across_Protocols.m first. To further analyze these results, one can also run Synergy_Across_Protocols_CompareFronts.m.
