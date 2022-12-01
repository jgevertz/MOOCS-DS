# MOOCS-DS
Multi-Objective Optimization of Combination Synergy - Dose Selection Algorithm
Data and code associated with "Guiding model-driven combination dose selection using multi-objective synergy optimization"
- All paper results related to the toy model are found by running MOOCS_DS_ToyModel.m in the ToyModel folder. Note you do need a version of MATLAB that can run a script for which the function files are embedded in the main m-file. If your version of MATLAB does not allow this, each function file needs to be put into its own .m file in the same directory and appropriately named.
- For the Pembro + Avastin model, run MOOCS_DS_PembroAvastin.m in the PembroAvastin folder. This creates all plots using the experimental protocol (Q3D). Note the same comment about embedded functions applies here. 
- To generate and visualize results for protocols other than the experimental protocol (for pembro + avastin model), run Synergy_Across_Protocols.m and/or Synergy_Across_Protocols_CompareFronts.m. Note you can run these without running MOOCS_DS_PembroAvastin.m, as the input files needed to make the visualizations are available in the CellLine1_H1299 and CellLine2_A549 folders. 
