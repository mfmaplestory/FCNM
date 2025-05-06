# Functional Connectivity Network Mapping (FCNM)

## Introduction
This repository contains MATLAB scripts and data for mapping functional connectivity networks. The project is designed to manage data from Excel to combined masks, individual-level FC maps, group-level t-maps, and final network probability maps. 

## Repository Structure
- `extra/`: Contains supplemental code and data necessary for running the main scripts.
- `publication_data/`: Contains partial intermediate and result data from the published articles.
- `FCNM.m`: Main script for processing data and generating network probability maps.
- `FCNM_Random_validation.m`: Script to validate that the network mappings are not random.
- `fc_analysis.m`: Parallel processing script for analyzing functional connectivity between ROIs and the whole brain.

## Prerequisites
- MATLAB (preferably R2019b or later)
- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
- Excels for input data preparation
- Your own time series data in 4D NIfTI format, as specified in the script (e.g., `Data_dir = 'F:\test\data\Hefei_GSR_656';`)

## Data Organization
Organize your data in the following structure:

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

**Important note:** The script expects each subject folder to contain a 4D BOLD time series file named specifically `swFiltered_4DVolume.nii`. If your files have different names, please modify the script accordingly.

## Running the Analysis
Run the main script:

```matlab
run fc_analysis.m
```

## Output Structure
The functional connectivity analysis toolbox generates the following outputs:

```
results/
└── fc_analysis/
    ├── roi_1/
    │   ├── roi_1.nii (copy of original ROI)
    │   ├── zFC_roi_1/ (Z-transformed connectivity maps)
    │   │   ├── zFC_subject_001.nii
    │   │   ├── zFC_subject_002.nii
    │   │   └── ...
    │   └── onesample_roi_1/ (Statistical results)
    │       ├── SPM.mat
    │       ├── spmT_0001.nii
    │       ├── spmT_0001_FWE_0.05.nii
    │       ├── spmT_0001_FWE_0.05_mask.nii
    │       ├── spmT_0001_positive_FDR_0.05.nii
    │       └── spmT_0001_positive_FDR_0.05_mask.nii
    └── ...
```

## Performance Tips
- Adjust the `num_workers` parameter to match your system's capabilities (default is 40)
- For very large datasets, consider using a cluster computing environment

## Publications: 
- Mo, Fan et al. "Network Localization of State and Trait of Auditory Verbal Hallucinations in Schizophrenia." Schizophrenia bulletin, sbae020. 24 Feb. 2024, doi:10.1093/schbul/sbae020
- Zhang, Xiaohan et al. "Brain Structural and Functional Damage Network Localization of Suicide." Biological psychiatry vol. 95,12 (2024): 1091-1099. doi:10.1016/j.biopsych.2024.01.003
- Cheng, Yan et al. "Brain network localization of gray matter atrophy, neurocognitive and social cognitive dysfunction in schizophrenia." Biological psychiatry, S0006-3223(24)01489-6. 3 Aug. 2024, doi:10.1016/j.biopsych.2024.07.021
