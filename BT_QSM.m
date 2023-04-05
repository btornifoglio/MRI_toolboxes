clear
clc
close all

te1 = .0046; %TE1 & TE spacing in [s] 
te = [te1:te1:6*te1]; % TE array

% load in MRIconvert .nii files of magnitude and phase
mag = load_nii('401_mag.nii'); 
mag1=mag.img;

phase = load_nii('401_phase.nii'); 
phase1=phase.img;

%% calculate r2star and echo-combined mag image (rms)
r2 = arlo(te,mag1);
echo_combined = rms(mag1,4);
niftiwrite(r2, '401_r2s.nii');
niftiwrite(echo_combined, '401_rms_echocombined.nii');
% use fslcpgeom, fslroi, and fslswapdim to get these to match rest
% %% saving .nii
% % trying to save into structure
% RtwoStar=mag;
% RtwoStar.img = r2; % cant save this properly so no point lol

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
[iFreq] = unwrapLaplacian(phase, matrix_size, voxel_size);
%% starting again
% after having exported dicoms from microDicom

metadata = dicoinfo("dicom-00001.dcm");
