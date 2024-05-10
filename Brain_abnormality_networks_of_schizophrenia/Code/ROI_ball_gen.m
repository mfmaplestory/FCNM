function ROI_ball_gen(coord,ROIs,ref_image,radius)
%% Format: ROI_ball_gen(coord,ROIs,ref_image,radius)
% =========================================================================
% generate binary ball masks according to the coordinations
% coord: vector of coordination which represents the center of the mask
% ROIs: cell list of output ROI name (.img or .nii postfix)
% ref_image: the reference image
% radius: the radius of the ball
%==========================================================================
clc

BAT_dir=which('BAT_fmri_batch');
BAT_dir=fileparts(BAT_dir);
switch nargin
    case 3
        radius=9;
    case 2
        radius=9;
        ref_image=[BAT_dir,filesep,'3mm_brainmask.nii'];
end
if nargin<2
    fprintf('Format: ROI_ball_gen(coord,ROIs,reference,radius)');
    fprintf('At least the coord and ROIs should be entered');
    exit
end

ref_vol=spm_vol(ref_image);
voxelsize=abs([ref_vol.mat(1,1),ref_vol.mat(2,2),ref_vol.mat(3,3)]);
origin=abs(ref_vol.mat(1:3,4)'./voxelsize);
ref_img=spm_read_vols(ref_vol);
m=size(ref_img);
[numROI, k0]=size(ROIs);
for r=1:numROI
    ROI_name=ROIs{r};
    mask=zeros(m);
    coord_mm=coord(r,:);
    coord_mm(1)=-coord_mm(1);
    coord_vox=round(coord_mm./voxelsize)+origin;
    radius_vox=round(radius/mean(voxelsize));
    xlim=coord_vox(1)-radius_vox:coord_vox(1)+radius_vox; fb=find(xlim>0&xlim<=m(1)); xlim=xlim(fb);
    ylim=coord_vox(2)-radius_vox:coord_vox(2)+radius_vox; fb=find(ylim>0&ylim<=m(2)); ylim=ylim(fb);
    zlim=coord_vox(3)-radius_vox:coord_vox(3)+radius_vox; fb=find(zlim>0&zlim<=m(3)); zlim=zlim(fb);
    for x=xlim
        for y=ylim
            for z=zlim
                Euclideand=pdist([coord_vox;[x y z]]);
                if Euclideand<=radius_vox
                    mask(x, y, z) =1;
                end
            end
        end
    end
    mvol=ref_vol;
    mvol.fname=[ROI_name,'.nii'];
    spm_write_vol(mvol,mask);
end
