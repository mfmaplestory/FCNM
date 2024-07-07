# Functional Connectivity Network Mapping (FCNM)

## Introduction
This repository contains MATLAB scripts and data for mapping functional connectivity networks. The project is designed to manage data from Excel to combined masks, individual-level FC maps, group-level t-maps, and final network probability maps. 

## Repository Structure
- `extra/`: Contains supplemental code and data necessary for running the main scripts.
- `publication_data/`: Contains partial intermediate and result data from the published articles.
- `FCNM.m`: Main script for processing data and generating network probability maps.
- `FCNM_Random_validation.m`: Script to validate that the network mappings are not random.

## Prerequisites
- MATLAB (preferably R2019b or later)
- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
- Excels for input data preparation
- Your own time series data in 4D NIfTI format, as specified in the script (e.g., `Data_dir = 'F:\test\data\Hefei_GSR_656';`)

## Publications: 
- Mo, Fan et al. “Network Localization of State and Trait of Auditory Verbal Hallucinations in Schizophrenia.” Schizophrenia bulletin, sbae020. 24 Feb. 2024, doi:10.1093/schbul/sbae020
- Zhang, Xiaohan et al. “Brain Structural and Functional Damage Network Localization of Suicide.” Biological psychiatry vol. 95,12 (2024): 1091-1099. doi:10.1016/j.biopsych.2024.01.003