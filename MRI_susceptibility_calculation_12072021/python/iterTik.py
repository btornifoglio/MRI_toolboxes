#DESCRIPTION: Calculates the susceptibility map by solving a least squares
#             problem with Tikhonov regularisation using the conjugate gradient method:
#             x = argmin_x || WM(B - B0(x*d)) ||_2^2 + alpha||Mx||_2^2 where * denotes
#             convolution
#
#Arguments:
#   field - The .nii file containing the field map (in ppm)
#   mask - The .nii file containing the binary tissue mask
#   noise - The .nii file containing the noise map (up to a scaling factor)
#   output - The .nii file where the susceptibility map is to be saved (in ppm)
#Optional arguments:
#   resolution - Image space resolution (default = 1,1,1)
#   B0direction - B0 direction (default = 0,0,1)
#   alpha - Regularisation parameter (default = 0.05)
#   stopping - Iteration stopping threshold (norm(res)/norm(b)) (default = 0.03 = 3%)
#
#AUTHOR:
#   Magnetic Resonance Imaging Group,
#   Department of Medical Physics and Biomedical Engineering,
#   University College London, UK, 2019

import os
import sys
import numpy as np
import nibabel as nib

def ConjugateGradient(A,b,tol):

    r = np.copy(b)
    p = np.copy(b)
    x = np.zeros(b.shape)

    while np.sqrt(np.dot(r,r)/np.dot(b,b))>tol:

        Ap = A(p)
        pAp = np.dot(p,Ap);
        alph = np.dot(p,r)/pAp

        x = x + alph*p
        r = r - alph*Ap

        Ar = A(r)
        bet = np.dot(p,Ar)/pAp
        p = r - bet*p

    return x

def SystemMatrix(kSpaceDipole,Mask,WeightsSquared,Alpha,x):
    x_shape = x.shape
    x = np.reshape(x,kSpaceDipole.shape)

    b = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))
    b = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.multiply(kSpaceDipole,b)))))

    b = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(np.multiply(WeightsSquared,b))))
    b = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.multiply(kSpaceDipole,b)))))

    b = b + Alpha*np.multiply(Mask,x)
    b = np.reshape(b,x_shape)

    return b

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
    print("Please specify the location of the binary tissue mask using the identifier: mask")
    sys.exit()

try:
    inputFile_noise = sys.argv[sys.argv.index("noise")+1]
    noise = nib.load(inputFile_noise)
    noise_data = noise.get_fdata()
except ValueError:
    print("Please specify the location of the noise map using the identifier: noise")
    sys.exit()

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
    Alpha = float(sys.argv[sys.argv.index("alpha")+1])
except ValueError:
    Alpha = 0.05

try:
    StoppingThreshold = float(sys.argv[sys.argv.index("stopping")+1])
except ValueError:
    StoppingThreshold = 0.03

# Calculate fitting Weights
Mask = np.copy(mask_data)

noise_data = np.multiply(noise_data,Mask)
ZerosInNoise = np.argwhere(noise_data==0)
noise_data[ZerosInNoise[:, 0], ZerosInNoise[:, 1], ZerosInNoise[:, 2]] = 1
Weights = np.true_divide(1,noise_data)
Weights[ZerosInNoise[:, 0], ZerosInNoise[:, 1], ZerosInNoise[:, 2]] = 0

Weights = Weights/(np.sum(Weights)/np.sum(Mask))
WeightsSquared = np.multiply(Weights,Weights)

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

# Calculate the right-hand side of the equations: b = DW^2 f
b = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(np.multiply(WeightsSquared,img_data))))
b = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.multiply(kSpaceDipole,b)))))
b = np.reshape(b,np.prod(b.shape))

# Define system matrix and solve using the conjugate gradient method
A = lambda x: SystemMatrix(kSpaceDipole,Mask,WeightsSquared,Alpha,x)
ImageSpaceSusc = ConjugateGradient(A,b,StoppingThreshold)
ImageSpaceSusc = np.reshape(ImageSpaceSusc,img_data.shape)

# PSF correction
kSpaceDipoleSquared = np.multiply(kSpaceDipole,kSpaceDipole)
PSF = np.real(np.fft.ifftn(np.fft.ifftshift(np.true_divide(kSpaceDipoleSquared,kSpaceDipoleSquared + Alpha))))
ImageSpaceSusc = np.multiply(ImageSpaceSusc/np.real(PSF[0,0,0]),Mask)

nib.save(nib.Nifti1Image(ImageSpaceSusc, np.eye(4)), outputFile)
