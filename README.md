# Influenza virus population cycles emerge from collections of variably responding cells
Jeremy D. Harris, Brigitte Martin, Molly E. Gallagher, Christopher B. Brooke, Katia Koelle July 2022

 -- *updated 07/11/22 by jdh* --

---

**Code for:** JD Harris, BE Martin, ME Gallagher, CB Brooke, KV Koelle. "Influenza virus population cycles emerge from collections of variably responding cells." 

- This repository contains all the MATLAB codes for simulating models and fitting models to data: <br>
  (1) Modelling and fitting to Genomic Equivalents (GE) for total viral yield in bulk cell culture;    <br>
  (2) The mean number of DIs equivalent to infection with at least 1 WT virus <br>
  (3) Modelling and fitting to Tissue Culture Infectious Dose 50\% (TCID$_{50}$) for the proportion of       WT virus in total viral yield.<br>
  (4) Simulating passage study models with best-fit parameter estimates

A preprint of the manuscript can be found on MedRxiv: [DOI](Not yet!) Not yet!!!

---

**Instructions:**
*MATLAB was used to parametrize models, run simulations, and plot figures. Once the Github repository is downloaded, navigate to the subdirectory 'Code_plt_ms_figures' to plot the figures in the manuscript by running the appropriate function. Simulation data can be reproduced by navigating to 'Code_sims' and running appropriate functions. See below for subfolder descriptions.*

---

**Folder descriptions:** <br>


- **Code_modelling:** All code for simulating model and fitting models to data: <br>
  (1) *'modelling_GE'*: fitting models of total viral yield vs. WT virus; estimating the mean number of DIs equivalent to infection with at least 1 WT virus <br>
  (2) *'modelling_passagestudy'*: simulate passage study models; create bifurcation diagram with cycle amplitude vs. DI fitness advantage  <br>
  (3) *'modelling_TCID50'*: fitting models of proportion WT in total viral yield to TCID50 vs. WT:DI MOI ratios <br>

Within each of these folders, the main script files parametrize and simulate the model. There are several choices for the user at the top of main file scripts. Typically, the first choice for the user is to save the simulation data using the variable 'save_ans': 0 means don't save and 1 means save. The output file 'filename' will be saved to the directory 'folder_location' so that the corresponding figure can be produced.

- **Code_plt_ms_figures:**
Each of the main scripts, reads in the data from 'Code_plt_ms_figures/data/' and the results  from 'Code_plt_ms_figures/results/' to plot the figures: Figures 1-7 and Supp Figures 1,2. If saved, the figures will be saved to the folder 'Figures_ms_all.' From here, they were uploaded to the Overleaf document.

- **Manuscript_forCodeReview:** Folder where you can find the draft of the manuscript that corresponds to Figures in 'Figures_ms_all.'

- **Figures_ms_all:** Figure files in the manuscript for main Figures 1-7 and Supplemental Figures 1,2.


- **Source_Data** Folder containing data shown and used: <br>
  - 'Source_Data_1.xlsx' - Passage study data: 2 lines over 73 passages; HAU and TCID50 measurements. <br>
  - 'Source_Data_2.xlsx' - Total viral yield (GE/mL) vs. WT mean MOI (WT stock) <br>
  - 'Source_Data_3.xlsx' - Total viral yield (GE/mL) and infectious virus (TCID$_{50}$) vs. WT:DI; the proportion WT virus by gene segment from WT stock/from DI-enriched stock  <br>

---

**References:** <br>
[1] <br>
[2] <br>
[3] <br>
[4]<br>
# cycling_dynamics
