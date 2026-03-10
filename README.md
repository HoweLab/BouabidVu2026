# Supplemental code for "An anatomical hotspot for striatal dopamine-acetylcholine interactions during reward and movement"

------------------
## Contents
* [Project Overview](#project-overview)
* [Data Organization](#data-organization)
* [Data Availability](#data-availability)
* [System Requirements](#system-requirements)
* [License and citation](#license-and-citation)
* [Acknowledgements](#acknowledgements)

------------------
## Project Overview

These scripts were custom-written by Mai-Anh Vu, except where indicated otherwise, to analyze the data in [Bouabid & Vu et al. (2026) An anatomical hotspot for striatal dopamine-acetylcholine interactions during reward and movement](https://www.biorxiv.org/content/10.64898/2026.01.20.700614v1).

The code included in this repository assumes data organization as detailed below in the Data Organization section. 
1. Preprocess:
   - batch_preprocess: calls preprocess_data to loop over data and run preprocess_data.
2. Analysis:
   - each of the scripts here can be run to execute the analyses relevant to the main figure indicated and will generate associated plots
   - functions called within each script are either included at the bottom of the script or in the common_functions folder
3. Analysis/common_functions:
   - functions shared across analyses
  

------------------
## Data Organization

The code included in this repository assume the following data organization:

1. Within the data folder are folders for each mouse, and within each mouse's folder are folders for the date of each experiment. In cases where multiple experiments were conducted on the same day, those are separated into folders with a letter subscript at the end. For example:
   - D:/UG27/240214
   - D:/UG27/240214r
   - D:/UG27/240220
     
2. Within each of the experiment folders are the extracted fluorescence files (has "ROIs" in the name) and behavior variable files (has "ttlIn1" or "ttlIn2" in the name). These data have been preprocessed (channels aligned, etc) into a single file with the file name format MOUSE_EXPDIR.mat, e.g., D:/UG27/240214r/UG27_240214r.mat. These files contains a struct with the following fields:
     - DA -- fluorescence data for DA
     - ACh -- fluorescence data for ACh
     - behav_DA -- behavior data aligned to DA fluorescence recording
     - behav_ACh -- behavior data aligned to ACh fluorescence recording
     - DA_idx -- the indices of the timepoints of the DA data that are included
     - ACh_idx -- the indices of the timepoints of the ACh data that are included
       
3. In the mouse's data folder is a spreadsheet fiber_table.xlsx that contains localization information for all of the recording locations

------------------
## Data Availability

Both raw datasets and processed and organized datasets are available as a [DANDISET](https://dandiarchive.org/dandiset/001692/0.260303.1724).

------------------
## System Requirements

Any operating system able of running Matlab R2020b. 

**Hardware requirement** : RAM (64GB), processor: Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz   3.19 GHz. 

**Software requirements**:  MATLAB can be installed on machines running Windows. MATLAB software requires a paid subscription to be used.

MATLAB2020b:
- This code was written and executed in MATLAB2020b. It has not been tested for compatibility with other versions. MATLAB software requires a paid subscription to be used.

MATLAB toolboxes:
- Curve Fitting Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
  
Other dependencies:
- [Red Blue Colormap](https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap)

------------------
##  License and Citation
 
If you use use this code in your research, please cite: 

Safa Bouabid*, Mai-Anh T. Vu*, Christian Noggle, Stefania Vietti-Michelina, Katherine Brimblecomb, Nicola Platt, Liangzhu Zhang, Anil Joshi, Stephanie Cragg, and Mark W. Howe.  (2026) [An anatomical hotspot for striatal dopamine-acetylcholine interactions during reward and movement](https://doi.org/10.64898/2026.01.20.700614).

This repository is released under the [MIT License](https://opensource.org/license/mit) - see the [LICENSE](LICENSE) file for details.

------------------
## Acknowledgements

This work was supported by the following funding sources: MWH - Aligning Science Across Parkinson’s **(ASAP-020370; ASAP-025192)** through the Michael J. Fox Foundation for Parkinson’s Research (MJFF), National Institute of Mental Health **(R01 MH125835)**, Whitehall Foundation Fellowship, Klingenstein-Simons Foundation Fellowship, Parkinson’s Foundation (Stanley Fahn Junior Faculty Award, **PF-SF-JFA-836662**); MTV - NIMH **F32MH120894**.

  
