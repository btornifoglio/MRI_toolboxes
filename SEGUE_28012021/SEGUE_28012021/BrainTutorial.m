
addpath('SEGUE_pfiles');

%% Use SEGUE for phase unwrapping

Inputs.Phase = your_phase_image; %3D array for single-echo, 4D array for multi-echo 
Inputs.Mask = your_tissue_mask; % 3D binary tissue mask, same size as one phase image

Unwrapped = SEGUE(Inputs);