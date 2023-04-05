
addpath('SEGUE_pfiles');

%% Use SEGUE separately in water and fat (or any other) masks

Inputs.Phase = load_nii_img_only('901_THRIVE_MP_TEs_1iso_PHASE.nii'); %3D array for single-echo, 4D array for multi-echo 
Inputs.Mask = load_nii_img_only('901_phaseunwrapMask.nii'); % e.g. water mask, 3D binary tissue mask, same size as one phase image
%Inputs.Mask(:,:,:,2) = second_tissue_mask; % e.g. fat mask, 3D binary tissue mask, same size as one phase image
% Inputs.Mask(:,:,:,3) = third_tissue_mask % etc.

Unwrapped = SEGUE(Inputs);

