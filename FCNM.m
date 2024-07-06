%%% ---------------------------------------------------------------------------------------------------
%%% Functional Connectivity Network Mapping (FCNM)
%%% ---------------------------------------------------------------------------------------------------
%%% This script manages data from Excel to combined masks, individual-level zFC,
%%% group-level T-maps, and final group probability maps.
%%% Developed by ChatGPT and Mofan, Anhui Medical University, December 28, 2022.
%%% ---------------------------------------------------------------------------------------------------
%%% Steps:
%%% 1. Load and set paths.
%%% 2. Establish directory structure and capture PDF file names for studies.
%%% 3. Prepare Excel spreadsheets for later ROI generation.
%%% 4. Generate ROI spheres, process BOLD files, and perform statistical analyses.
%%% 5. Create and refine group probability maps.
%%% ---------------------------------------------------------------------------------------------------

%% Step 1: Add spm12 and ROI_ball_gen_combined.m to the path and clear all variables

clear, clc
support_path = 'F:\test\extra'; % Path for support files
addpath(support_path);

%% Step 2: Read the basic information and title of each independent study
% Typically, we have downloaded all the PDFs of included studies, recommended to be named as "P + Number + Author + Year of Publication", e.g., P01Li2017
% If not all PDFs are downloaded, create a placeholder PDF for missing studies as we only need to read each document's basic information and title

path = 'F:\test';
mkdir([path,'\Articles_Included']); % Create a new folder 'Articles_Included' to store all included PDF files of the studies
Articles_path = [path, '\Articles_Included'];
File = dir(fullfile(Articles_path, '*.pdf')); % Read all PDF files in 'Articles_Included', representing all studies. The file name without '.pdf' is obtained by (end-4)
Filename = {File.name};


%% Step 3: After reading the information, create an Excel spreadsheet for all articles' coordinates, which will be manually filled in later

mkdir([path, filesep, 'Articles_Included_Excel']); % Create a new folder 'Articles_Included_Excel'
cd([path, filesep, 'Articles_Included_Excel']); % Enter this folder

for i = 1:length(Filename) % Create an Excel spreadsheet for each study, where coordinates will be manually entered into the spreadsheet later and saved. If there are multiple ROI results, the spreadsheet can be autofilled downwards
    filename = char(Filename(i));
    xlswrite([filename(1:end-4), '.xlsx'], cellstr([filename(1:end-4), '_ROI01']), 'sheet1', 'D1')
end

%% Step 4: Generate spheres for each study's ROI and merge spheres from the same study

radius_Seeds_results = [path, '\4mm_Seeds_Results\']; % Note to adjust the radius size when generating spheres with different radii, here it is 4mm

parfor i = 1:length(Filename)
    filename = char(Filename(i));
    article_folder = fullfile(radius_Seeds_results, filename(1:end-4));
    mkdir(article_folder)

    ROI_Seeds_path = fullfile(article_folder, 'ROI_Seeds');
    if ~exist(ROI_Seeds_path, 'dir')
        mkdir(ROI_Seeds_path)
    end

    onesample_path = fullfile(article_folder, ['onesample_', filename(1:end-4)]);
    if ~exist(onesample_path, 'dir')
        mkdir(onesample_path)
    end

    zFC_path = fullfile(article_folder, ['zFC_', filename(1:end-4)]);
    if ~exist(zFC_path, 'dir')
        mkdir(zFC_path)
    end

    % Proceed to generate and combine ROI spheres
    [coord, ROIs] = xlsread([path, filesep, 'Articles_Included_Excel', filesep, filename(1:end-4), '.xlsx']);
    output_path = ROI_Seeds_path;
    combined_mask = ROI_ball_gen_combined(coord, ROIs, [support_path, filesep, 'BNA_mask_3m.nii'], 4, output_path);


    % Loop through each participant's directory to process BOLD files
    Data_dir = 'F:\test\data\Hefei_GSR_656';  % Specify the path where the 4D NIfTI time series data are stored
    participant_folders = dir(fullfile(Data_dir, '*'));
    participant_folders = participant_folders(~ismember({participant_folders.name}, {'.', '..'}));

    Mask_BNA = spm_read_vols(spm_vol(fullfile(support_path, 'BNA_mask_3m.nii')));
    Mask_Brain = spm_read_vols(spm_vol(fullfile(support_path, 'BrainMask_05_61x73x61.img')));

    for p = 1:length(participant_folders)

        boldFiles = dir(fullfile(Data_dir, participant_folders(p).name, '*.nii'));
        boldFile = fullfile(boldFiles.folder, boldFiles.name);
        bold_vol = spm_vol(boldFile);
        boldData = spm_read_vols(bold_vol);
        boldData = boldData .* Mask_Brain;
        boldData2D = reshape(boldData, [], size(boldData, 4));


        roiTimeSeries = mean(boldData2D(combined_mask > 0, :), 1);
        maskIndex = Mask_Brain > 0;
        boldDataMasked = boldData2D(maskIndex(:), :);

        fcMap = corr(roiTimeSeries', boldDataMasked');
        zfcMap = zeros(size(boldData, 1), size(boldData, 2), size(boldData, 3));
        zfcMap(maskIndex) = atanh(fcMap);


        zFC_filename = fullfile(zFC_path, ['zFC_', participant_folders(p).name, '.nii']);
        zFC_vol = bold_vol(1);
        zFC_vol.fname = zFC_filename;
        spm_write_vol(zFC_vol, zfcMap);

    end

    zFC_files = spm_select('FPList', zFC_path, '^zFC_.*\.nii$');

    if isempty(zFC_files)
        warning('No zFC files found in %s.', zFC_path);
        continue;
    end

    % Setup the statistical design for one-sample t-tests in SPM
    matlabbatch = [];
    matlabbatch{1}.spm.stats.factorial_design.dir = {onesample_path};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(zFC_files);
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(support_path, 'BNA_mask_3m.nii,1')};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    % Model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(onesample_path, 'SPM.mat')};

    % Contrast specification
    matlabbatch{3}.spm.stats.con.spmmat = {fullfile(onesample_path, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'One Sample T-test';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    % Run the batch job
    spm_jobman('run', matlabbatch);


    % Load SPM T-map file
    V = spm_vol(fullfile(onesample_path, 'spmT_0001.nii'));
    [X, XYZ] = spm_read_vols(V);

    % Compute p-values for positive T statistics
    P_uncorrected = 1 - spm_Tcdf(X(X > 0), 655);
    P_sorted = sort(P_uncorrected);
    P_corrected = spm_P_FDR(P_sorted);

    % Reorder corrected P-values
    [~, sort_idx] = sort(P_uncorrected, 'ascend');
    [~, unsort_idx] = sort(sort_idx, 'ascend');
    P_corrected = P_corrected(unsort_idx);

    % Apply FDR threshold
    thresholded_map = zeros(size(X));
    positive_voxel_indices = find(X > 0);
    thresholded_map(positive_voxel_indices(P_corrected < 0.05)) = X(positive_voxel_indices(P_corrected < 0.05));

    % Save the corrected T-map
    V.fname = fullfile(onesample_path, 'spmT_0001_positive_FDR.nii');
    spm_write_vol(V, thresholded_map);

    % Create a new mask image for voxels above threshold
    mask = thresholded_map > 0;
    V.fname = fullfile(onesample_path, 'spmT_0001_positive_FDR_mask.nii');
    spm_write_vol(V, mask);

end

%% Step 5: Create and refine group probability maps.

output_folder_avg = fullfile(radius_Seeds_results, 'Average_probability_map');
mkdir(output_folder_avg);

% Initialize the accumulation matrix
sum_matrix = [];

% Retrieve paths to all 'onesample_*' subfolders
onesample_folders = dir(fullfile(radius_Seeds_results, '*', 'onesample_*'));
num_folders = length(onesample_folders);

% Iterate through each 'onesample_*' folder and accumulate the contents of each mask file
for i = 1:num_folders
    parent_folder_name = onesample_folders(i).folder;
    onesample_folder_name = onesample_folders(i).name;
    mask_file = dir(fullfile(parent_folder_name, onesample_folder_name, '*_positive_FDR_mask.nii'));

    if ~isempty(mask_file)
        file_path = fullfile(parent_folder_name, onesample_folder_name, mask_file.name);
        finalVol = spm_vol(file_path);
        [Y, ~] = spm_read_vols(finalVol);

        if isempty(sum_matrix)
            sum_matrix = zeros(size(Y));
        end

        sum_matrix = sum_matrix + Y;
    end
end

% Calculate the average
avg_matrix = sum_matrix / num_folders;

% Retain values greater than or equal to 0.6
avg_matrix(avg_matrix < 0.6) = 0;

% Remove clusters smaller than 30 voxels
% Use the bwlabeln function to label connected regions
[L, num] = bwlabeln(avg_matrix, 6);

% Iterate through each connected region and remove those smaller than 30 voxels
for j = 1:num
    if sum(L(:) == j) < 30
        avg_matrix(L == j) = 0;
    end
end

% Save the final averaged matrix to a new nii file
final_mask_filename = fullfile(output_folder_avg, 'average_mask_60_30.nii');
finalVol.fname = final_mask_filename;
spm_write_vol(finalVol, avg_matrix);