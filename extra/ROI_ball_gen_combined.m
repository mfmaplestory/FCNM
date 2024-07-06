function combined_mask = ROI_ball_gen_combined(coord, ROIs, ref_image, radius, output_path)
%% Format: ROI_ball_gen(coord, ROIs, ref_image, radius, output_path)
% =========================================================================
% generate binary ball masks according to the coordinations
% coord: vector of coordination which represents the center of the mask
% ROIs: cell list of output ROI name (.img or .nii postfix)
% ref_image: the reference image
% radius: the radius of the ball
% output_path: the path where the combined seeds file will be saved
% output: combined balls mask, representing merged ROIs.
% =========================================================================
clc

if nargin < 5  % Check if the output path is provided
    output_path = ''; % Default to current directory or specify a default path
end

BAT_dir = which('BAT_fmri_batch');
BAT_dir = fileparts(BAT_dir);
switch nargin
    case 3
        radius = 9;
    case 2
        radius = 9;
        ref_image = [BAT_dir, filesep, '3mm_brainmask.nii'];
end
if nargin < 2
    fprintf('Format: ROI_ball_gen(coord, ROIs, reference, radius, output_path)');
    fprintf('At least the coord and ROIs should be entered');
    exit
end

ref_vol = spm_vol(ref_image);
voxelsize = abs([ref_vol.mat(1,1), ref_vol.mat(2,2), ref_vol.mat(3,3)]);
origin = abs(ref_vol.mat(1:3,4)' ./ voxelsize);
ref_img = spm_read_vols(ref_vol);
m = size(ref_img);
combined_mask = zeros(m);  % Initialize the total mask
[numROI, ~] = size(ROIs);

for r = 1:numROI
    mask = zeros(m);
    coord_mm = coord(r,:);
    coord_mm(1) = -coord_mm(1);
    coord_vox = round(coord_mm ./ voxelsize) + origin;
    radius_vox = round(radius / mean(voxelsize));
    xlim = coord_vox(1) - radius_vox : coord_vox(1) + radius_vox;
    fb = find(xlim > 0 & xlim <= m(1)); xlim = xlim(fb);
    ylim = coord_vox(2) - radius_vox : coord_vox(2) + radius_vox;
    fb = find(ylim > 0 & ylim <= m(2)); ylim = ylim(fb);
    zlim = coord_vox(3) - radius_vox : coord_vox(3) + radius_vox;
    fb = find(zlim > 0 & zlim <= m(3)); zlim = zlim(fb);
    
    for x = xlim
        for y = ylim
            for z = zlim
                Euclideand = pdist([coord_vox; [x y z]]);
                if Euclideand <= radius_vox
                    mask(x, y, z) = 1;
                end
            end
        end
    end
    combined_mask = combined_mask + mask;  % Add the current mask to the total mask
end

combined_mask(combined_mask > 1) = 1;  % Set all values greater than 1 to 1

% Save the combined mask
mvol = ref_vol;
output_filename = fullfile(output_path, [ROIs{1}(1:end-6), '_Seeds.nii']);  % Using fullfile for proper path construction
mvol.fname = output_filename;
spm_write_vol(mvol, combined_mask);

end

