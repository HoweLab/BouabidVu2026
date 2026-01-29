# BouabidVu2026
Analysis code to accompany Bouabid, Vu, et al., 2026: https://www.biorxiv.org/content/10.64898/2026.01.20.700614v1

## Contents
The code included in this repository assume data organization as detailed below
1. preprocess:
   - batch_preprocess: calls preprocess_data to loop over data and preprocess
2. analysis:
   - each of the scripts here can be run to execute the analyses relevant to the main figure indicated and will generate associated plots
   - functions called within each script are either included at the bottom of the script, or, for functions shared across scripts, included in the common_functions folder
3. analysis/common_functions:
   - functions shared across analyses

## Data Organization:
The code included in this repository assume the following data organization:
1. Within the data folder are folders for each mouse, and within each mouse's folder are folders for the date of each experiment. In cases where multiple experiments were conducted on the same day, those are separated into folders with a letter subscript at the end. For example:
- D:/UG27/240214
- D:/UG27/240214r
- D:/UG27/240220
2. Within each of the experiment folders are the extracted fluorescence files (has "ROIs" in the name) and behavior variable files (has "ttlIn1" or "ttlIn2" in the name). Since these are dual DA-ACh recordings, these data have been aligned into a single file with the file name format MOUSE_EXPDIR.mat, e.g., D:/UG27/240214r/UG27_240214r.mat. These files contains a struct with the following fields:
  - DA -- fluorescence data for DA
  - ACh -- fluorescence data for ACh
  - behav_DA -- behavior data aligned to DA fluorescence recording
  - behav_ACh -- behavior data aligned to ACh fluorescence recording
  - DA_idx -- the indices of the timepoints of the DA data that are included
  - ACh_idx -- the indices of the timepoints of the ACh data that are included
3. In the mouse's data folder is a spreadsheet fiber_table.xlsx that contains localization information for all of the recording locations
