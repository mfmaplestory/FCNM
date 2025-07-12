# Functional Connectivity Network Mapping (FCNM)

## Introduction
This repository contains MATLAB scripts and data for mapping functional connectivity networks. The project is designed to manage data from Excel to combined masks, individual-level FC maps, group-level t-maps, and final network probability maps.

## Repository Structure
* `extra/`: Contains supplemental code and data necessary for running the main scripts.
* `publication_data/`: Contains partial intermediate and result data from the published articles.
* `FCNM.m`: Main script for processing data and generating network probability maps from Excel coordinates.
* `FCNM_Random_validation.m`: Script to validate that the network mappings are not random.
* `ROI_FC_fwe__fdr_parfor.m`: Supplementary script for optimized parallel processing of existing ROIs for whole-brain functional connectivity analysis.

## Prerequisites
* MATLAB (preferably R2019b or later)
* SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
* Excels for input data preparation (for FCNM.m)
* Your own time series data in 4D NIfTI format, as specified in the script (e.g., `Data_dir = 'F:\test\data\Hefei_GSR_656';`)

## Usage Notes
The repository contains two complementary analysis approaches:
1. **FCNM.m** - Use this script when starting from Excel coordinates to generate ROIs and then perform connectivity analysis.
2. **ROI_FC_fwe__fdr_parfor.m** - Use this script when you already have ROI masks and want to perform optimized parallel functional connectivity analysis.

## Enhanced Functional Connectivity Analysis

For users with existing ROI masks, the `ROI_FC_fwe__fdr_parfor.m` script offers optimized parallel processing for functional connectivity analysis. This script can significantly improve processing speed when analyzing multiple ROIs.

### Data Organization for ROI_FC_fwe__fdr_parfor.m
```
data/
├── roi_masks/
│   ├── roi_1.nii
│   ├── roi_2.nii
│   └── ...
└── bold_data/
    ├── subject_001/
    │   └── swFiltered_4DVolume.nii  # 4D BOLD time series for subject 1
    ├── subject_002/
    │   └── swFiltered_4DVolume.nii  # 4D BOLD time series for subject 2
    └── ...
```

**Note:** The script expects each subject folder to contain a 4D BOLD time series file named `swFiltered_4DVolume.nii`.

### Output Structure
```
results/
└── fc_analysis/
    ├── roi_1/
    │   ├── roi_1.nii (copy of original ROI)
    │   ├── zFC_roi_1/ (Z-transformed connectivity maps)
    │   └── onesample_roi_1/ (Statistical results)
    │       ├── SPM.mat
    │       ├── spmT_0001.nii
    │       ├── spmT_0001_FWE_0.05.nii
    │       ├── spmT_0001_FWE_0.05_mask.nii
    │       ├── spmT_0001_positive_FDR_0.05.nii
    │       └── spmT_0001_positive_FDR_0.05_mask.nii
    └── ...
```

## Publications
* Mo, Fan et al. "Network Localization of State and Trait of Auditory Verbal Hallucinations in Schizophrenia." Schizophrenia bulletin, sbae020. 24 Feb. 2024, doi:10.1093/schbul/sbae020

* Cheng, Yan et al. "Brain network localization of gray matter atrophy, neurocognitive and social cognitive dysfunction in schizophrenia." Biological psychiatry, S0006-3223(24)01489-6. 3 Aug. 2024, doi:10.1016/j.biopsych.2024.07.021

* Yao, Shanwen et al. "Network localization of genetic risk for schizophrenia and bipolar disorder."
