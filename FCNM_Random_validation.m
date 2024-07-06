%%%----------------------------------------------------------------------------------------------------
%%%---------------------------------------------------------------------------------------------------
%%%---Supplement for FCNM----------------------------------------------------------------------------
%%%---Validating Random Networks---1000 iterations---40 cores---Extended runtime--------------------
%%%---Output is 1000 nii images---yielding 1000 final result images---probability maps-------------
%%%---ChatGPT---Mofan---Anhui Medical University---2023.12.24---------------------------------------
%%%---BNA_mask_3m---Excel_example---See extra file---------------------------------------------------
%%%---------------------------------------------------------------------------------------------------
%%%----------------------------------------------------------------------------------------------------

clear, clc

% Set the number of cores for the parallel pool
numCores = 40;

% Start the parallel pool at the beginning of the code
pool = gcp('nocreate'); % Check if a parallel pool already exists
if isempty(pool)
    parpool(numCores); % Start a parallel pool with the specified number of cores
end

% Set the number of iterations
num_iterations = 1000;

% Create the main output folder
boldPath = 'D:\J_Team\Grade4\mofan\Hefei_GSR';
total_output_folder = 'D:\J_Team\Grade4\mofan\Random_Validation\Total_Results_calzFC';
mkdir(total_output_folder);

% Load NIfTI file (Mask) and Excel file---see the example attachment on GitHub
[Mask_BNA, ~, ~, Header] = y_ReadAll('D:\J_Team\Grade4\mofan\Random_Validation\BNA_mask_3m.nii');
[~, txt, raw] = xlsread('D:\J_Team\Grade4\mofan\Random_Validation\Excel_example.xlsx');

% Retrieve paths of all BOLD files
boldDirs = dir(fullfile(boldPath, '**', '*.nii'));
boldFiles = fullfile({boldDirs.folder}, {boldDirs.name});

% Begin loop
parfor iteration = 1: num_iterations
    % Create a unique directory for each iteration
    iter_folder = fullfile(total_output_folder, sprintf('Iteration_%d', iteration));
    mkdir(iter_folder);

    % ------ Part 1: Process ROIs and create masks ------

    % Create a directory to save Seeds masks
    output_folder_seeds = fullfile(iter_folder, 'Seeds_mask');
    mkdir(output_folder_seeds);

    % Process each ROI
    for i = 1:size(raw, 1)
        % Name and number of points for the ROI
        roi_name = raw{i, 1};
        num_points = raw{i, 2};

        % Find coordinates of all points with a value of 1
        [x, y, z] = ind2sub(size(Mask_BNA), find(Mask_BNA > 0));

        % Initialize a matrix of zeros
        new_ROI = zeros(size(Mask_BNA));

        % Randomly select the specified number of points
        rand_indices = randperm(length(x), num_points);
        selected_points = [x(rand_indices), y(rand_indices), z(rand_indices)];

        % Set the value of each selected point and its immediate neighbors to 1
        for j = 1:num_points
            px = selected_points(j, 1);
            py = selected_points(j, 2);
            pz = selected_points(j, 3);

            % Set the value to 1 for the selected point and its direct neighbors
            for dx = -1:1
                for dy = -1:1
                    for dz = -1:1
                        % Ensure inclusion of only direct neighbors (front, back, left, right, up, down)
                        if abs(dx) + abs(dy) + abs(dz) == 1 || (dx == 0 && dy == 0 && dz == 0)
                            if px+dx > 0 && px+dx <= size(Mask_BNA, 1) && ...
                                    py+dy > 0 && py+dy <= size(Mask_BNA, 2) && ...
                                    pz+dz > 0 && pz+dz <= size(Mask_BNA, 3)
                                new_ROI(px+dx, py+dy, pz+dz) = 1;
                            end
                        end
                    end
                end
            end
        end

        % Save the new ROI matrix to the specified path
        out_name = fullfile(output_folder_seeds, [roi_name '.nii']);
        y_Write(new_ROI, Header, out_name);

        % ------ Part 2: Process each Mask and run statistical analysis ------

        % Path to the Seeds mask folder
        First_path = output_folder_seeds;
        Second_path = fullfile(iter_folder, 'Seeds_Results');
        mkdir(Second_path);

        % Process each Mask
        File = dir(fullfile(First_path, '*.nii'));
        SingleFilename = File(i).name(1:end-4);
        Path_basic = fullfile(Second_path, SingleFilename);
        mkdir(Path_basic);
        copyfile(fullfile(First_path, [SingleFilename '.nii']), Path_basic);

        % Create required sub-folders
        mkdir(fullfile(Path_basic, ['onesample_', SingleFilename]));
        mkdir(fullfile(Path_basic, ['zFC_', SingleFilename]));


        % Loop to process each BOLD file
        for k = 1:length(boldFiles)
            % Read BOLD data and apply mask
            [boldData, ~, ~, boldHeader] = y_ReadAll(boldFiles{k});
            boldData = boldData .* Mask_BNA;

            % Reshape boldData into a 2D array
            boldData2D = reshape(boldData, [], size(boldData, 4));

            % Calculate the average time series for the ROI
            roiTimeSeries = mean(boldData2D(new_ROI > 0, :), 1);

            % Select voxels of interest using maskData
            maskIndex = Mask_BNA > 0;
            boldDataMasked = boldData2D(maskIndex(:), :);

            % Calculate FC map (vectorized operation)
            fcMap = corr(roiTimeSeries', boldDataMasked');

            % Reshape the FC map to match original BOLD data dimensions
            zfcMap = zeros(size(boldData, 1), size(boldData, 2), size(boldData, 3));
            zfcMap(maskIndex) = atanh(fcMap); % Apply Fisher z transformation

            % Save zFC map
            boldHeader.fname = fullfile(Path_basic, ['zFC_', SingleFilename], ['zFC_', num2str(k), '.nii']);
            y_Write(zfcMap, boldHeader, boldHeader.fname);
        end


        % Code to generate SPM statistical results
        zFC_namelist = dir(fullfile(Path_basic, ['zFC_', SingleFilename], '*.nii'));
        len = length(zFC_namelist);
        file_name = cell(1, len);
        for k = 1:len
            file_name{k} = fullfile(Path_basic, ['zFC_', SingleFilename], zFC_namelist(k).name);
        end
        File_name = file_name(3:end)';


        jobs = struct();

        % Set parameters for factorial_design
        jobs.stats{1}.spm.stats.factorial_design.dir = {fullfile(Path_basic, ['onesample_', SingleFilename])};
        jobs.stats{1}.spm.stats.factorial_design.des.t1.scans = File_name;
        jobs.stats{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        jobs.stats{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        jobs.stats{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        jobs.stats{1}.spm.stats.factorial_design.masking.im = 1;
        jobs.stats{1}.spm.stats.factorial_design.masking.em = {'D:\J_Team\Grade4\mofan\Random_Validation\BNA_mask_3m.nii,1'};
        jobs.stats{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        jobs.stats{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        jobs.stats{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        % Set parameters for model estimation
        jobs.stats{2}.spm.stats.fmri_est.spmmat = {fullfile(Path_basic, ['onesample_', SingleFilename], 'SPM.mat')};
        jobs.stats{2}.spm.stats.fmri_est.write_residuals = 0;
        jobs.stats{2}.spm.stats.fmri_est.method.Classical = 1;

        % Set parameters for contrast specification
        jobs.stats{3}.spm.stats.con.spmmat = {fullfile(Path_basic, ['onesample_', SingleFilename], 'SPM.mat')};
        jobs.stats{3}.spm.stats.con.consess{1}.tcon.name = 'onesample';
        jobs.stats{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        jobs.stats{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        jobs.stats{3}.spm.stats.con.delete = 0;

        % Run SPM jobs
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs.stats);


        rmdir(fullfile(Path_basic, ['zFC_', SingleFilename]), 's');
        % Read SPM T-map files
        V = spm_vol(fullfile(Path_basic, ['onesample_', SingleFilename], 'spmT_0001.nii'));
        [X, XYZ] = spm_read_vols(V);

        % Get P-values for positive T-statistics
        P_uncorrected = 1 - spm_Tcdf(X(X > 0), 655);

        % Sort P-values
        P_sorted = sort(P_uncorrected);

        % Calculate FDR-corrected P-values
        P_corrected = spm_P_FDR(P_sorted);

        % To match corrected P-values to their corresponding voxels, indexing is needed
        [~, sort_idx] = sort(P_uncorrected, 'ascend');
        [~, unsort_idx] = sort(sort_idx, 'ascend');
        P_corrected = P_corrected(unsort_idx);

        % Apply FDR < 0.05 threshold
        thresholded_map = zeros(size(X));
        positive_voxel_indices = find(X > 0);
        thresholded_map(positive_voxel_indices(P_corrected < 0.05)) = X(positive_voxel_indices(P_corrected < 0.05));

        % Save the corrected T-map
        V.fname = fullfile(Path_basic, ['onesample_', SingleFilename], 'spmT_0001_positive_FDR.nii');
        spm_write_vol(V, thresholded_map);

        % Create a new mask image, setting voxels greater than 0 to 1
        mask = thresholded_map > 0;
        V.fname = fullfile(Path_basic, ['onesample_', SingleFilename],  'spmT_0001_positive_FDR_mask.nii');
        spm_write_vol(V, mask);
        %toc
    end

    % ------ Part 3: Accumulate and average Masks, remove small voxel clusters ------

    output_folder_avg = fullfile(iter_folder, 'Average_50_30');
    mkdir(output_folder_avg);

    % Initialize the accumulation matrix
    sum_matrix = [];

    % Retrieve paths to all 'onesample_*' subfolders
    onesample_folders = dir(fullfile(Second_path, '*', 'onesample_*'));
    num_folders = length(onesample_folders);

    % Iterate through each 'onesample_*' folder and accumulate the contents of each mask file
    for i = 1:num_folders
        parent_folder_name = onesample_folders(i).folder;
        onesample_folder_name = onesample_folders(i).name;
        mask_file = dir(fullfile(parent_folder_name, onesample_folder_name, '*_positive_FDR_mask.nii'));

        if ~isempty(mask_file)
            file_path = fullfile(parent_folder_name, onesample_folder_name, mask_file(1).name);
            V = spm_vol(file_path);
            [Y, ~] = spm_read_vols(V);

            if isempty(sum_matrix)
                sum_matrix = zeros(size(Y));
            end

            sum_matrix = sum_matrix + Y;
        end
    end

    % Calculate the average
    avg_matrix = sum_matrix / num_folders;

    % Retain values greater than or equal to 0.5
    avg_matrix(avg_matrix < 0.5) = 0;

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
    V.fname = fullfile(output_folder_avg,[SingleFilename, '_average_mask_50_30.nii']);
    spm_write_vol(V, avg_matrix);

    % Reset variables at the end of iteration
    avg_matrix = []; % Reset avg_matrix to an empty array
    V = struct();    % Reset V to a new empty structure
end
