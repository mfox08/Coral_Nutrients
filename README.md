# Differential Acclimation to chronic nutrient enrichment 

This project examines the physilogical effects of chronic nutrient encirhment on 2 species of reef-building coral. The corals Pocillopora acuta and Porites compressa were in Kāne'ohe Bay, O'ahu, Hawai'i and acclimated to 5 levels of nutrient enrichment for a 5-week period.

Folders:
Data
Manuscript
Code

Data - this folder contains all the raw data files required by the code files and will recreate all analyses and figures for the manuscript. 
* Isotope nubbins.csv - the list of all corals analyzed for stable isotopes.
* Isotopes_raw.csv - raw stable isotope data that is merged with metadata from the above file.
* Kaneohe_historical_nutrients.csv - baseline nutrient concentrations for the study site from Drupp et al. 2011.
* Kbay_nutrients_experimental.csv - weekly measurements of nutrients and water chemistry from experimental aquaria.
* PAM_master_file.csv - master file of all maximum quantum yield measurements for the experiment.
* Phys_master_file.csv - master file with all physiological parameters measured during the study.
* Tank_irradiance.csv - ambient irradiance data recorded in the center of the experimental tank set up over the experiment.
* Tank_temperatures.csv - water temperature of experimental tanks throughout the experiment.

    * This group of files is used for the analysis of the respirometry data in the physiology code file
        * Buoyant_weights_11.24.2015.csv - bouyant weights of all coral fragments used to correct incumation chamber water volumes
        * CHAIN_respirometry_chamberID_mid_final.csv - metadata for each incubation chamber 
        * CHAIN_Slopes.csv - raw oxygen production data
        * Nubbin_growth_and_SA_05.02.2016.csv - surface area measurements for all coral framgents for standardization
        * TPAIN_master_all_variables09.17.csv - metadata for all corals used for the respirometry measurements 


manuscript
PDF of online supporting information

Code:
1. Fox_etal2020_tank_conditions.R - analysis of environmental conditions during the course of the experiment
2. Fox_etal2020_nutrient_analyses.R - analysis of historical baseline and experimental nutrient concentrations
3. Fox_etal2020_Physiology_analyses.R - analysis of all physiological and photosynthetic data and creation of main text figs. 1,2,3
4. Fox_etal2020_Isotope_analyses.R - analysis of stable isotope data and creation of main text Fig. 4

