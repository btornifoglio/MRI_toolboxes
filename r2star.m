clear
clc
close all

te1 = .0046; %te and te spacing; in seconds 
te = [te1 2*te1 3*te1 4*te1 5*te1 6*te1];

mag = load_nii('901_THRIVE_MP_TEs_1iso_MAG.nii'); 
mag1=mag.img;

phase = load_nii('901_THRIVE_MP_TEs_1iso_PHASE.nii'); 
phase1=phase.img;

%% calculate r2star and echo-combined mag image (rms)
r2 = arlo(te,mag1);
echo_combined = rms(mag1,4);


%% saving .nii
% resaving all as coordinate system is different from fsleyes
niftiwrite(r2,'r2star.nii');
niftiwrite(echo_combined,'echoCombined.nii');
niftiwrite(mag1, 'mag.nii');
niftiwrite(phase1, 'phase.nii')

%% echo to echo inconsistencies

% Remove echo-to-echo phase inconsistencies in readout phase corected
% complex data
% NOTE: initial testing suggests that use of this function DOES NOT
% introduce negative effects further down the pipeline if normal data
% (i.e., without phase correction) is processed with this function.

%[iField] = iField_correction(iField,voxel_size);
% 
%% phase unwrap

matrix_size = size(phase1); 
matrix_size = matrix_size(1:2)
voxel_size = [0.525 0.525 0.5];
[iFreq ] = unwrapLaplacian(phase, matrix_size, voxel_size);