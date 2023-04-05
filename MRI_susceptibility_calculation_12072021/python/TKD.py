#DESCRIPTION: Calculates the susceptibility map from the input field map based on
#             Ref. Shmueli et al. MRM 62, 2009 and Ref. Schweser et al. MRM, 2012
#
#Arguments:
#   field - The .nii file containing the field map (in ppm)
#   output - The .nii file where the susceptibility map is to be saved (in ppm)
#Optional arguments:
#   mask - The .nii file containing the binary tissue mask (default = 1)
#   resolution - Image space resolution (default = 1,1,1)
#   B0direction - B0 direction (default = 0,0,1)
#   threshold - k-space truncation threshold (default = 2/3)
#
#AUTHOR:
#   Magnetic Resonance Imaging Group,
#   Department of Medical Physics and Biomedical Engineering,
#   University College London, UK, 2019

import os
import sys
import numpy as np
import nibabel as nib

# Sort input parameters
try:
    inputFile_field = sys.argv[sys.argv.index("field")+1]
    img = nib.load(inputFile_field)
    img_data = img.get_fdata()
    img_hdr = img.header
except ValueError:
    print("Please specify the location of the input field map using the identifier: field")
    sys.exit()

try:
    inputFile_mask = sys.argv[sys.argv.index("mask")+1]
    mask = nib.load(inputFile_mask)
    mask_data = mask.get_fdata()
except ValueError:
    mask_data = np.ones(img_data.shape)

try:
    outputFile = sys.argv[sys.argv.index("output")+1]
except ValueError:
    print("Please specify the location of the output susceptibility map using the identifier: output")
    sys.exit()

try:
    resolution = np.fromstring(sys.argv[sys.argv.index("resolution")+1], dtype=int, sep=',')
except ValueError:
    resolution = [1,1,1]

try:
    B0direction = np.fromstring(sys.argv[sys.argv.index("B0direction")+1], dtype=int, sep=',')
except ValueError:
    B0direction = [0,0,1]

try:
    Threshold = float(sys.argv[sys.argv.index("threshold")+1])
except ValueError:
    Threshold = 2/3

# Transform data into Fourier space
kSpaceimg = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(img_data)))

# Create dipole in k-space
matrixSize = np.array(img_data.shape)
kSpaceResolution = np.divide(1,np.multiply(resolution,matrixSize))
kSpaceEdge1 = -np.floor(matrixSize/2)
kSpaceEdge2 = kSpaceEdge1 + matrixSize

Y, X, Z = np.meshgrid(np.arange(kSpaceEdge1[1],kSpaceEdge2[1])*kSpaceResolution[1],
 np.arange(kSpaceEdge1[0],kSpaceEdge2[0])*kSpaceResolution[0],
 np.arange(kSpaceEdge1[2],kSpaceEdge2[2])*kSpaceResolution[2])

kSquared = np.multiply(X,X) + np.multiply(Y,Y) + np.multiply(Z,Z)
ZerosInkSquared = np.argwhere(kSquared==0)
kSquared[ZerosInkSquared[:, 0], ZerosInkSquared[:, 1], ZerosInkSquared[:, 2]] = 1

kSpaceDipole = 1/3 - np.true_divide(np.square(B0direction[0]*X +
B0direction[1]*Y + B0direction[2]*Z),kSquared)
kSpaceDipole[ZerosInkSquared[:, 0], ZerosInkSquared[:, 1], ZerosInkSquared[:, 2]] = 0

# Threshold kernel
InverseKernel = np.copy(kSpaceDipole)
ZerosInkSpaceDipole = np.argwhere(InverseKernel==0)
InverseKernel[ZerosInkSpaceDipole[:, 0], ZerosInkSpaceDipole[:, 1], ZerosInkSpaceDipole[:, 2]] = np.inf

LessThanThreshold = np.argwhere(np.absolute(InverseKernel)<Threshold)
InverseKernel[LessThanThreshold[:, 0], LessThanThreshold[:, 1], LessThanThreshold[:, 2]] = Threshold*\
np.sign(InverseKernel[LessThanThreshold[:, 0], LessThanThreshold[:, 1], LessThanThreshold[:, 2]])

# Apply kernel
kSpaceSusc = np.true_divide(kSpaceimg,InverseKernel)
ImageSpaceSusc = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(kSpaceSusc))))

# PSF correction
PSF = np.real(np.fft.ifftn(np.fft.ifftshift(np.true_divide(kSpaceDipole,InverseKernel))))
ImageSpaceSusc = np.multiply(ImageSpaceSusc/np.real(PSF[0,0,0]),mask_data)

nib.save(nib.Nifti1Image(ImageSpaceSusc, np.eye(4)), outputFile)
