%%%----------------------------------------------------------------------------------------------------
%%%---------------------------------------------------------------------------------------------------
%%%---Functional_connectivity_network_mapping---FCNM----------------------------------------------------------------------------------
%%%---The output result is the group-level functional connectivity T-map--------------------------------------------
%%%---Mofan---Anhui_Medical_University---2022.10.28----------------------------------------
%%%---------------------------------------------------------------------------------------------------
%%%----------------------------------------------------------------------------------------------------

%% Step 1: Add spm12, dpabi, and ROI_ball_gen.m (used for drawing spheres) to the path and clear all variables
% Add spm12, dpabi, and ROI_ball_gen.m manually
clear, clc
addpath('C:\Matlab\Matlab_zuoye\coordinate_ROI'); %Add the path for ROI_ball_gen.m

%% Step 2: Read the basic information and title of each independent study
% Typically, we have downloaded all the PDFs of included studies, recommended to be named as "P + Number + Author + Year of Publication", e.g., P01Li2017
% If not all PDFs are downloaded, create a placeholder PDF for missing studies as we only need to read each document's basic information and title

path = 'F:\test'; % Modify the path as needed, this is the main folder path for our project, e.g., 'Disease_Network_Mapping'. Adjust according to your own file structure
mkdir([path, filesep, 'Articles_Included']); % Create a new folder 'Articles_Included' to store all included PDF files of the studies
Articles_path = [path, filesep, 'Articles_Included'];
File = dir(fullfile(Articles_path, '*.pdf')); % Read all PDF files in 'Articles_Included', representing all studies. The file name without '.pdf' is obtained by (end-4)
Filename = {File.name};
[~, col] = size(Filename);
%% Step 3: After reading the information, create an Excel spreadsheet for all articles' coordinates, which will be manually filled in later

mkdir([path, filesep, 'Articles_Included_Excel']); % Create a new folder 'Articles_Included_Excel'
cd([path, filesep, 'Articles_Included_Excel']); % Enter this folder

for i = 1:col % Create an Excel spreadsheet for each study, where coordinates will be manually entered into the spreadsheet later and saved. If there are multiple ROI results, the spreadsheet can be autofilled downwards
    filename = char(Filename(i));
    xlswrite([filename(1:end-4), '.xlsx'], cellstr([filename(1:end-4), '_ROI01']), 'sheet1', 'D1')
end

%% Step 4: Generate spheres for each study's ROI and merge spheres from the same study

[row, col] = size(Filename);
for i = 1:col
    filename = char(Filename(i));
    radius_Seeds_Excel = [path, '\4mm_Seeds_Excel\']; % Note to adjust the radius size when generating spheres with different radii, here it is 4mm
    mkdir([radius_Seeds_Excel, filename(1:end-4)]) % Create a folder for each article
    cd([radius_Seeds_Excel, filename(1:end-4)]);
    mkdir ROI_Seeds
    mkdir(['onesample_', filename(1:end-4)])
    mkdir(['zFC_', filename(1:end-4)])
    % Draw individual ROI spheres
    [coord, ROIs] = xlsread([path, filesep, 'Articles_Included_Excel', filesep, filename(1:end-4), '.xlsx']); % Read coordinates from the Excel file
    Reference_image_path = 'F:\test'; % Modify the path as needed, reference file required for drawing ROI spheres should match the resolution of all subjects' images
    ROI_ball_gen(coord, ROIs, [Reference_image_path, filesep, 'BNA_mask_3m.nii'], 4); % 'BNA_mask_3m.nii' is the reference image for ROI resolution; 4 is the sphere radius of 4mm;

    % Merge individual ROI spheres into a large mask for each study
    Combine_file = {};
    for j = 1:length(ROIs)
        combine_file = [radius_Seeds_Excel, filename(1:end-4), '\', filename(1:end-4), '_ROI', num2str(j, '%02d'), '.nii,1'];
        Combine_file = [Combine_file; {combine_file}];
    end
    Combine_ROI = [];
    for k = 1:length(ROIs)
        combine_ROI = ['i', int2str(k)];
        Combine_ROI = [Combine_ROI, '+', combine_ROI];
        Combine_ROIs = Combine_ROI(2:end);
    end
    
    jobs{1}.spm.util.imcalc.input = Combine_file;
    jobs{1}.spm.util.imcalc.output = [filename(1:end-4), '_Seeds'];
    jobs{1}.spm.util.imcalc.outdir = {[radius_Seeds_Excel, filename(1:end-4), '\ROI_Seeds']};
    jobs{1}.spm.util.imcalc.expression = Combine_ROIs;
    jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    jobs{1}.spm.util.imcalc.options.dmtx = 0;
    jobs{1}.spm.util.imcalc.options.mask = 0;
    jobs{1}.spm.util.imcalc.options.interp = 1;
    jobs{1}.spm.util.imcalc.options.dtype = 4;
    spm('defaults', 'FMRI')
    spm_jobman('run', jobs)
    clear jobs
    % Use spm to merge Masks and rename
end

%% Step 5: Begin functional connectivity analysis with DPABI
working_directory = 'W:\'; % Define the main folder location for all subjects
% Note, we need to manually load a parameter file from DPABI software. Only for the first step, the rest is automated
load('F:\test\dpabi.mat'); % Locate and load this mat file
for i = 1:col
    filename = char(Filename(i));
    Cfg.CalFC.ROIDef{1, 1} = [radius_Seeds_Excel, filename(1:end-4), '\ROI_Seeds\', filename(1:end-4), '_Seeds.nii'];
    
    AutoDataProcessParameter = Cfg;
    WorkingDir = [];
    SubjectListFile = [];
    IsAllowGUI = [];
    [Error, AutoDataProcessParameter] = DPARSFA_run(AutoDataProcessParameter, WorkingDir, SubjectListFile, IsAllowGUI); % This step begins the calculation of functional connectivity
    
    file = dir([working_directory, filesep, 'Results\', 'FC_', Cfg.StartingDirName, filesep, 'zFC*']);
    
    for j = 1:length(file)
        movefile([working_directory, filesep, 'Results\', 'FC_', Cfg.StartingDirName, filesep, file(j).name], ...
        [radius_Seeds_Excel, filename(1:end-4), '\zFC_', filename(1:end-4)]);
    end
    rmdir([working_directory, filesep, 'Masks'], 's');
    rmdir([working_directory, filesep, 'Results'], 's'); % After copying the results, they can be deleted
    movefile([radius_Seeds_Excel, filename(1:end-4), '\P*'], [radius_Seeds_Excel, filename(1:end-4), '\ROI_Seeds']) % Organize files

% Now, start the one-sample T-test
% Till_spmT_0001
    filename = char(Filename(i));
    zFC_namelist = dir([radius_Seeds_Excel, filename(1:end-4), '\zFC_', filename(1:end-4), '\zFC*']);
    len = length(zFC_namelist);
    for j = 1:len
        file_name{j} = [radius_Seeds_Excel, filename(1:end-4), '\zFC_', filename(1:end-4), '\', zFC_namelist(j).name];
    end
    File_name = file_name(1:end)';
  %-----------------------------------------------------------------------
  % Factorial design specification
  jobs{1}.spm.stats.factorial_design.dir = {[radius_Seeds_Excel, filename(1:end-4), '\onesample_', filename(1:end-4)]};
  jobs{1}.spm.stats.factorial_design.des.t1.scans = File_name;
  jobs{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
  jobs{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  jobs{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  jobs{1}.spm.stats.factorial_design.masking.im = 1;
  jobs{1}.spm.stats.factorial_design.masking.em = {'F:\test\BNA_mask_3m.nii,1'}; % Note, modify to the path of the reference mask
  jobs{1}.spm.stats.factorial_design.globalc.g_omit = 1;
  jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
  jobs{1}.spm.stats.factorial_design.globalm.glonorm = 1;
  
  spm('defaults', 'FMRI');
  spm_jobman('run', jobs);
  clear jobs
  
  
  jobs{1}.spm.stats.fmri_est.spmmat = {[radius_Seeds_Excel, filename(1:end-4), '\onesample_', filename(1:end-4), '\SPM.mat']};
  jobs{1}.spm.stats.fmri_est.write_residuals = 0;
  jobs{1}.spm.stats.fmri_est.method.Classical = 1;
  
  spm('defaults', 'FMRI');
  spm_jobman('run', jobs);
  clear jobs
  % Model estimation
  
  jobs{1}.spm.stats.con.spmmat = {[radius_Seeds_Excel, filename(1:end-4), '\onesample_', filename(1:end-4), '\SPM.mat']};
  jobs{1}.spm.stats.con.consess{1}.tcon.name = 'onesample';
  jobs{1}.spm.stats.con.consess{1}.tcon.weights = 1;
  jobs{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
  jobs{1}.spm.stats.con.delete = 0;
  
  spm('defaults', 'FMRI');
  spm_jobman('run', jobs);
  clear jobs
end

% For subsequent processing of the generated T-maps, further steps such as correction, binarization, overlaying, and averaging can be customized according to individual requirements