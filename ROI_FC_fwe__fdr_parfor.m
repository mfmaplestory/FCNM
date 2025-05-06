%% Functional Connectivity Analysis Script (Parallel Version)
% This script processes multiple ROI masks and BOLD data, calculates functional connectivity,
% and performs statistical analysis using SPM.
%
% Author: Äª·² (Mo Fan)
% Version: 1.0.0
% GitHub: https://github.com/mfmaplestory/FCNM
%
% Features:
% - Parallel processing of multiple ROIs for improved efficiency
% - Fisher's Z transformation of correlation maps
% - Statistical analysis with one-sample t-tests
% - FWE and FDR correction methods
%
% Requirements:
% - MATLAB (tested on R2019b or later)
% - SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
%
% Citation:
% No citation needed - this code was co-created by human and AI partners.

%% Clear workspace
clear;
clc;

%% CONFIG - Edit these parameters according to your needs
% =========================================================================
% Path settings - CHANGE THESE PATHS to your local directories
% -------------------------------------------------------------------------
config = struct();

% ROI mask folder path - folder containing your ROI mask files (.nii)
config.roi_folder = './data/roi_masks'; 

% BOLD data main folder path
% This folder should contain subject subfolders, where each subject folder contains
% a 4D BOLD time series file named 'swFiltered_4DVolume.nii'
% Structure: bold_main_folder/subject_001/swFiltered_4DVolume.nii
config.bold_main_folder = './data/bold_data';

% Output results folder path
config.output_folder = './results/fc_analysis';

% Brain mask path (for analysis restriction)
config.brain_mask_path = './template/BrainMask_05_61x73x61.nii';

% Parallel computing settings
% -------------------------------------------------------------------------
config.use_parallel = true;    % Whether to use parallel computing
config.num_workers = min(40, feature('numcores')); % Adjust based on CPU cores

% Statistical analysis settings
% -------------------------------------------------------------------------
config.fwe_threshold = 0.05;   % FWE correction threshold
config.fdr_threshold = 0.05;   % FDR correction threshold

% Display settings
% -------------------------------------------------------------------------
config.verbose = true;         % Display detailed processing information
config.log_interval = 10;      % Log progress every N subjects

%% Create output directories
% Create main output folder if it doesn't exist
if ~exist(config.output_folder, 'dir')
    mkdir(config.output_folder);
    if config.verbose
        fprintf('Created output directory: %s\n', config.output_folder);
    end
end

%% Initialize parallel computing pool
if config.use_parallel
    % Check if parallel pool already exists, if not create one
    if isempty(gcp('nocreate'))
        pool = parpool('local', config.num_workers);
        if config.verbose
            fprintf('Created parallel pool with %d workers\n', pool.NumWorkers);
        end
    else
        pool = gcp('nocreate');
        if config.verbose
            fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
        end
    end
end

%% Get subject folders and ROI files
% Get all subject folders
bold_subfolders = get_subject_folders(config.bold_main_folder);
if config.verbose
    fprintf('Found %d subject folders\n', length(bold_subfolders));
end

% Get all ROI mask files
roi_files = dir(fullfile(config.roi_folder, '*.nii'));
if config.verbose
    fprintf('Found %d ROI mask files\n', length(roi_files));
end

%% Load brain mask
brain_mask = load_brain_mask(config);
if config.verbose
    fprintf('Brain mask loaded or created\n');
end

%% Create necessary folders for each ROI
for roi_idx = 1:length(roi_files)
    create_roi_folders(roi_idx, roi_files, config);
end

%% Process each ROI mask (parallel or serial)
if config.use_parallel
    if config.verbose
        fprintf('Processing ROIs using parallel computing...\n');
    end
    parfor roi_idx = 1:length(roi_files)
        % Direct function call for parallel execution
        process_single_roi(roi_idx, roi_files, config, bold_subfolders, brain_mask);
    end
else
    if config.verbose
        fprintf('Processing ROIs sequentially...\n');
    end
    for roi_idx = 1:length(roi_files)
        process_single_roi(roi_idx, roi_files, config, bold_subfolders, brain_mask);
    end
end

%% Clean up parallel pool
if config.use_parallel && ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
    if config.verbose
        fprintf('Closed parallel computing pool\n');
    end
end

fprintf('All ROI functional connectivity analyses completed!\n');

%% Helper Functions

% Get subject folders
function subfolders = get_subject_folders(main_folder)
    subfolders = dir(main_folder);
    subfolders = subfolders([subfolders.isdir]); % Keep only directories
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Remove . and ..
end

% Load or create brain mask
function brain_mask = load_brain_mask(config)
    if exist(config.brain_mask_path, 'file')
        % Use SPM to read the brain mask file
        V = spm_vol(config.brain_mask_path);
        brain_mask = spm_read_vols(V);
    else
        % If brain mask not provided, create a temporary one from first BOLD data
        fprintf('Brain mask not found, creating a temporary mask from first BOLD data\n');
        first_subj = dir(fullfile(config.bold_main_folder, '*'));
        first_subj = first_subj([first_subj.isdir]);
        first_subj = first_subj(~ismember({first_subj.name}, {'.', '..'}));
        
        if ~isempty(first_subj)
            first_bold_path = fullfile(config.bold_main_folder, first_subj(1).name, 'swFiltered_4DVolume.nii');
            if exist(first_bold_path, 'file')
                % Read the first volume of the 4D BOLD data
                V = spm_vol([first_bold_path, ',1']); % Only read the first volume
                temp_bold = spm_read_vols(V);
                brain_mask = ones(size(temp_bold));
            else
                error('Cannot find BOLD data for temporary brain mask creation');
            end
        else
            error('No subject folders found to create a temporary brain mask');
        end
    end
end

% Create folders for ROI analysis
function create_roi_folders(roi_idx, roi_files, config)
    roi_name = roi_files(roi_idx).name(1:end-4); % Remove .nii extension
    
    % Create current ROI output folder
    roi_output_folder = fullfile(config.output_folder, roi_name);
    if ~exist(roi_output_folder, 'dir')
        mkdir(roi_output_folder);
    end
    
    % Create zFC subfolder
    zfc_folder = fullfile(roi_output_folder, ['zFC_', roi_name]);
    if ~exist(zfc_folder, 'dir')
        mkdir(zfc_folder);
    end
    
    % Create statistical analysis subfolder
    stats_folder = fullfile(roi_output_folder, ['onesample_', roi_name]);
    if ~exist(stats_folder, 'dir')
        mkdir(stats_folder);
    end
    
    % Copy ROI mask to output folder
    roi_path = fullfile(config.roi_folder, roi_files(roi_idx).name);
    copyfile(roi_path, roi_output_folder);
end

% Process a single ROI
function zfc_file_list = process_single_roi(roi_idx, roi_files, config, bold_subfolders, brain_mask)
    roi_name = roi_files(roi_idx).name(1:end-4); % Remove .nii extension
    if config.verbose
        fprintf('Processing ROI: %s (%d/%d)\n', roi_name, roi_idx, length(roi_files));
    end
    
    % Get current ROI output folder
    roi_output_folder = fullfile(config.output_folder, roi_name);
    zfc_folder = fullfile(roi_output_folder, ['zFC_', roi_name]);
    
    % Read current ROI mask using SPM
    roi_path = fullfile(config.roi_folder, roi_files(roi_idx).name);
    V_roi = spm_vol(roi_path);
    roi_mask = spm_read_vols(V_roi);
    
    % Save all subject zFC image paths for later statistical analysis
    zfc_file_list = {};
    
    %% Process each subject
    for subj_idx = 1:length(bold_subfolders)
        subj_name = bold_subfolders(subj_idx).name;
        
        % Log progress at intervals
        if config.verbose && (mod(subj_idx, config.log_interval) == 1 || subj_idx == length(bold_subfolders))
            fprintf('  ROI: %s - Processing subject progress: %d/%d\n', roi_name, subj_idx, length(bold_subfolders));
        end
        
        % Get current subject's BOLD data path
        bold_path = fullfile(config.bold_main_folder, subj_name, 'swFiltered_4DVolume.nii');
        
        % Skip if BOLD file doesn't exist
        if ~exist(bold_path, 'file')
            warning('BOLD data not found for subject %s', subj_name);
            continue;
        end
        
        % Read BOLD data using SPM
        V_bold = spm_vol(bold_path);
        bold_data = spm_read_vols(V_bold);
        
        % Calculate functional connectivity map
        zfc_volume = calculate_fc(bold_data, roi_mask, brain_mask);
        
        % Save zFC image using SPM
        zfc_file_name = ['zFC_', subj_name, '.nii'];
        zfc_file_path = fullfile(zfc_folder, zfc_file_name);
        
        % Create a new header for the zFC volume
        V_out = V_bold(1); % Copy header from first volume
        V_out.fname = zfc_file_path;
        V_out.dt = [16, 0]; % Set to float32 data type
        V_out.n = [1 1]; % Only one volume
        
        % Write the zFC volume
        spm_write_vol(V_out, zfc_volume);
        
        % Add current zFC file path to list
        zfc_file_list{end+1} = zfc_file_path;
    end
    
    % Perform statistical analysis
    perform_statistical_analysis(zfc_file_list, roi_name, roi_output_folder, config);
end

% Calculate functional connectivity
function zfc_volume = calculate_fc(bold_data, roi_mask, brain_mask)
    % Get data dimensions
    [nx, ny, nz, nt] = size(bold_data);
    
    % Reshape 4D BOLD data to 2D matrix (voxels ¡Á time points)
    bold_data_2d = reshape(bold_data, [], nt);
    
    % Calculate ROI's mean time series
    roi_indices = find(roi_mask > 0);
    roi_time_series = mean(bold_data_2d(roi_indices, :), 1);
    
    % Select voxels within brain mask
    brain_indices = find(brain_mask > 0);
    brain_time_series = bold_data_2d(brain_indices, :);
    
    % Calculate correlation between ROI and each brain voxel
    fc_map = zeros(length(brain_indices), 1);
    
    % Optimized correlation calculation - vectorized approach
    % Normalize the ROI time series
    roi_ts_norm = (roi_time_series - mean(roi_time_series)) / std(roi_time_series);
    
    % Normalize each voxel's time series
    brain_ts_norm = zeros(size(brain_time_series));
    for i = 1:size(brain_time_series, 1)
        brain_ts_norm(i,:) = (brain_time_series(i,:) - mean(brain_time_series(i,:))) / std(brain_time_series(i,:));
    end
    
    % Calculate correlation using dot product of normalized time series
    for i = 1:length(brain_indices)
        fc_map(i) = (roi_ts_norm * brain_ts_norm(i,:)') / (nt - 1);
    end
    
    % Apply Fisher's Z transformation
    zfc_map = atanh(fc_map);
    
    % Reshape FC map to 3D volume
    zfc_volume = zeros(nx, ny, nz);
    zfc_volume(brain_indices) = zfc_map;
end

% Perform statistical analysis
function perform_statistical_analysis(zfc_file_list, roi_name, roi_output_folder, config)
    % Statistical analysis subfolder path
    stats_folder = fullfile(roi_output_folder, ['onesample_', roi_name]);
    
    %% Perform statistical analysis (one-sample t-test)
    if length(zfc_file_list) > 1
        if config.verbose
            fprintf('Performing statistical analysis for ROI %s\n', roi_name);
        end
        
        % Setup SPM batch job
        jobs = setup_spm_jobs(zfc_file_list, stats_folder, config);
        
        % Run SPM batch job
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs.stats);
        
        % Perform FWE correction
        perform_fwe_correction(stats_folder, zfc_file_list, roi_name, config);
        
        % Perform FDR correction
        perform_fdr_correction(stats_folder, zfc_file_list, roi_name, config);
    else
        warning('Not enough zFC images for statistical analysis for ROI %s', roi_name);
    end
end

% Setup SPM batch jobs
function jobs = setup_spm_jobs(zfc_file_list, stats_folder, config)
    jobs = struct();
    
    % Factorial design (one-sample t-test)
    jobs.stats{1}.spm.stats.factorial_design.dir = {stats_folder};
    jobs.stats{1}.spm.stats.factorial_design.des.t1.scans = zfc_file_list';
    jobs.stats{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    jobs.stats{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    jobs.stats{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    jobs.stats{1}.spm.stats.factorial_design.masking.im = 1;
    jobs.stats{1}.spm.stats.factorial_design.masking.em = {config.brain_mask_path};
    jobs.stats{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    jobs.stats{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    jobs.stats{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Model estimation
    jobs.stats{2}.spm.stats.fmri_est.spmmat = {fullfile(stats_folder, 'SPM.mat')};
    jobs.stats{2}.spm.stats.fmri_est.write_residuals = 0;
    jobs.stats{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % Contrast setup
    jobs.stats{3}.spm.stats.con.spmmat = {fullfile(stats_folder, 'SPM.mat')};
    jobs.stats{3}.spm.stats.con.consess{1}.tcon.name = 'onesample';
    jobs.stats{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    jobs.stats{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    jobs.stats{3}.spm.stats.con.delete = 0;
    
end

% Perform FWE correction
function perform_fwe_correction(stats_folder, zfc_file_list, roi_name, config)
    if config.verbose
        fprintf('  Performing FWE correction (p < %.3f) for ROI %s...\n', config.fwe_threshold, roi_name);
    end
    
    % Read T map
    V = spm_vol(fullfile(stats_folder, 'spmT_0001.nii'));
    [stat_map, ~] = spm_read_vols(V);
    
    % Load SPM.mat to get necessary parameters
    load(fullfile(stats_folder, 'SPM.mat'));
    
    % Calculate FWE-corrected threshold
    df = SPM.xX.erdf; % Residual degrees of freedom
    R = SPM.xVol.R;   % RESEL counts
    S = sum(SPM.xVol.S(:)); % Search volume
    
    % Calculate whole-brain FWE-corrected threshold (p < config.fwe_threshold)
    intensity_threshold = spm_uc(config.fwe_threshold, [1, df], 'T', R, 1, S);
    
    % Apply threshold
    fwe_map = zeros(size(stat_map));
    fwe_map(stat_map > intensity_threshold) = stat_map(stat_map > intensity_threshold);
    
    % Create mask (binary image)
    mask_map = zeros(size(stat_map));
    mask_map(fwe_map > 0) = 1;
    
    % Save results
    V.fname = fullfile(stats_folder, sprintf('spmT_0001_FWE_%.3f.nii', config.fwe_threshold));
    spm_write_vol(V, fwe_map);
    
    V.fname = fullfile(stats_folder, sprintf('spmT_0001_FWE_%.3f_mask.nii', config.fwe_threshold));
    spm_write_vol(V, mask_map);
    
    % Output results information
    if config.verbose
        if any(mask_map(:))
            fprintf('  ROI %s FWE correction completed (p < %.3f), threshold T = %.3f\n', ...
                roi_name, config.fwe_threshold, intensity_threshold);
            fprintf('  Found %d significant voxels\n', sum(mask_map(:)));
        else
            fprintf('  ROI %s FWE correction completed, but no significant voxels found (p < %.3f, FWE)\n', ...
                roi_name, config.fwe_threshold);
        end
    end
end

% Perform FDR correction
function perform_fdr_correction(stats_folder, zfc_file_list, roi_name, config)
    % Read SPM T-map
    t_map_path = fullfile(stats_folder, 'spmT_0001.nii');
    V = spm_vol(t_map_path);
    [X, ~] = spm_read_vols(V);
    
    % Get uncorrected p-values for positive values
    df = length(zfc_file_list) - 1; % Degrees of freedom
    P_uncorrected = 1 - spm_Tcdf(X(X > 0), df);
    
    % Sort P-values
    P_sorted = sort(P_uncorrected);
    
    % Calculate FDR-corrected P-values
    P_corrected = spm_P_FDR(P_sorted);
    
    % Match corrected P-values to corresponding voxels
    [~, sort_idx] = sort(P_uncorrected, 'ascend');
    [~, unsort_idx] = sort(sort_idx, 'ascend');
    P_corrected = P_corrected(unsort_idx);
    
    % Apply FDR < config.fdr_threshold threshold
    thresholded_map = zeros(size(X));
    positive_voxel_indices = find(X > 0);
    sig_indices = positive_voxel_indices(P_corrected < config.fdr_threshold);
    
    if ~isempty(sig_indices)
        thresholded_map(sig_indices) = X(sig_indices);
        
        % Save corrected T-map
        V.fname = fullfile(stats_folder, sprintf('spmT_0001_positive_FDR_%.3f.nii', config.fdr_threshold));
        spm_write_vol(V, thresholded_map);
        
        % Create binary mask image
        mask = thresholded_map > 0;
        V.fname = fullfile(stats_folder, sprintf('spmT_0001_positive_FDR_%.3f_mask.nii', config.fdr_threshold));
        spm_write_vol(V, mask);
        
        if config.verbose
            fprintf('  ROI %s statistical analysis completed, found %d significant voxels (FDR < %.3f)\n', ...
                roi_name, sum(mask(:)), config.fdr_threshold);
        end
    else
        if config.verbose
            fprintf('  ROI %s statistical analysis completed, no significant voxels found (FDR < %.3f)\n', ...
                roi_name, config.fdr_threshold);
        end
    end
end