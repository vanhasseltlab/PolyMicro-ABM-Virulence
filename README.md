# PolyMicro-ABM-Virulence
## An agent based model of a polymicrobial community with virulence-mediated interspecies interactions

This folder contains all scripts and processed data needed to reproduce the analysis and figures shown in the manuscript “Investigating the dynamics of species-specific virulence traits in polymicrobial infections”. All code was produced by Eske van Meegen and Catharina Herzberg, Leiden Academic Center for Drug Research, Leiden University, The Netherlands. Parts of the model (Eske_Model_v6.py) and data processing is based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo.  doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM.

## Requirements:
The analysis has been performed with python 3.8. The python package mesa 1.1.1 was used to create the model. Further dependencies: matplotlib 3.5.3, numpy 1.23.3, pandas 1.4.4

## Overview of files: 

The analysis was performed in 3 steps: 1. simulation, 2. data processing 3. plotting. The following lists the files in the order of execution.

## Simulation scripts:

    1. Virulence strengths simulations scripts (02-03):
        ◦ Input: Eske_Model_v5.py, Eske_PostSim_v4.py
    2. Virulence simulations scripts (04-10):
        ◦ Input: Eske_Model_v6.py, Eske_PostSim_v4.py
    3. Antibiotic simulation scripts (20-22)
        ◦ Input: Eske_Model_v6.py, Eske_PostSim_v4.py

## Processing scrips:

    1. Create complete database files from each set of simulations:
        ◦ Create complete database “00_Complete_Database”: Eske_Database_make_complete_database.py
        ◦ Create complete database for virulence strengths data “001_Complete_Database_Virulence_Strengths”: Eske_Database_make_complete_database_VirulenceStrengths.py
    2. Calculate AUC and add to complete databases
        ◦ for “00_Complete_Database”: Eske_Database_calculate_auc.py
        ◦ for “001_Complete_Database_Virulence_Strengths”: Eske_Database_calculate_auc_VirulenceStrengths.py
          
## Plotting scripts:

    1. Plot Figures 3, 4 & 5: Plot1_Virulence_dynamics.py
        ◦ from Data“001_Complete_Database_Virulence_Strengths“
        ◦ figures saved in “Figures/VirulenceDynamics”
    2. Plot Figures 6 & 7: Plot2_AUC_Maps.py 
        ◦ from Data “00_Complete_Database”
        ◦ figures saved in “Figures/AUCs_together”
    3. Plot Figures 8 & 9: Plot3_Growthcurves.py (also filles up database pre-plotting) 
        ◦ from “00_Complete_Database”
        ◦ figures saved in “Figures/Growthcurves”
