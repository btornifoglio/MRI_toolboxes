close all
clear
clc

%DicomFolder = '/Volumes/tornifob/3T-Philips/230301-Ph/rawdata/901-classicD/';
%DicomList = dir('*.dcm');

%[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,files]=Read_DICOM_Philips_3D(DicomFolder);

wrappedPhase = load_nii_img_only('901_THRIVE_MP_TEs_1iso_PHASE.nii');
mask = load_nii_img_only('901_phaseunwrapMask.nii');
matrixSize = [400 400];
voxelSize = [0.525 0.525];

unwrappedPhase = UnwrapPhaseMacro(wrappedPhase,matrixSize,voxelSize,'method','laplacian');