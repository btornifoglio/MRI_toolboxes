function SusceptibilityMap = TKD(Parameters)
%DESCRIPTION: SusceptibilityMap = TKD(Parameters)
%             Calculates the susceptibility map from the input field map based on 
%             Ref. Shmueli et al. MRM 62, 2009 and Ref. Schweser et al. MRM, 2012 
%
%INPUTS:
%   Parameters(struct): Parameters.FieldMap(double array) - input image in ppm
%                       Parameters.Mask(array) - binary tissue mask (default = 1)
%                       Parameters.Threshold(double) - kernel threshold (default = 2/3)
%                       Parameters.Resolution(double vector) - image resolution vector (dx,dy,dz) in mm (default = 1 mm isotropic)
%                       Parameters.B0direction(double vector) - 3-element unit vector aligned with B0 (default = [0,0,1]) 
%
%OUTPUTS:
%   SusceptibilityMap(double array) - output susceptibility map in ppm
%
%DEPENDENCIES:
%   None
%
%AUTHOR: 
%   Magnetic Resonance Imaging Group, 
%   Department of Medical Physics and Biomedical Engineering, 
%   University College London, UK, 2019

% Sort input parameters
if isfield(Parameters,'FieldMap')
    FieldMap = double(Parameters.FieldMap);
else
    warndlg('Please specify Parameters.FieldMap!','Warning')
    return;
end
if isfield(Parameters,'Mask')
    Mask = double(Parameters.Mask);
else
    Mask = ones(size(FieldMap));
end
if isfield(Parameters,'Threshold')
    Threshold = Parameters.Threshold;
else
    Threshold = 2/3;
end
if isfield(Parameters,'Resolution')
    Resolution = Parameters.Resolution;
else
    Resolution = [1 1 1];
end
if isfield(Parameters,'B0direction')
    B0direction = Parameters.B0direction;
else
    B0direction = [0 0 1];
end

% Transform data into Fourier space
kSpaceField = fftshift(fftn(ifftshift(FieldMap)));

% Create dipole in k-space
MatrixSize = size(FieldMap);
dkx = 1/Resolution(1)/MatrixSize(1);
dky = 1/Resolution(2)/MatrixSize(2);
dkz = 1/Resolution(3)/MatrixSize(3);

kSpaceEdge1 = -floor(MatrixSize/2);
kSpaceEdge2 = kSpaceEdge1 + MatrixSize - 1;

[Y,X,Z] = meshgrid((kSpaceEdge1(2):kSpaceEdge2(2))*dky,...
                   (kSpaceEdge1(1):kSpaceEdge2(1))*dkx,...
                   (kSpaceEdge1(3):kSpaceEdge2(3))*dkz);

kSquared = X.^2+Y.^2+Z.^2;
               
kSpaceDipole = 1/3 - (B0direction(1)*X + B0direction(2)*Y + B0direction(3)*Z).^2./kSquared;
kSpaceDipole(isnan(kSpaceDipole)) = 0;

% Create kernel
Kernel = 1./kSpaceDipole;
LessThanThreshold = find(abs(kSpaceDipole)<Threshold);
Kernel(LessThanThreshold) = (1/Threshold)*sign(kSpaceDipole(LessThanThreshold));
Kernel(kSpaceDipole==0) = 0;

% Calculate PSF correction factor
PSF = ifftn(ifftshift(kSpaceDipole.*Kernel));

% 'Deconvolution'
kSpaceSusceptibility = kSpaceField.*Kernel/(PSF(1));

% Inverse Fourier transform
SusceptibilityMap = fftshift(ifftn(ifftshift(kSpaceSusceptibility)));
SusceptibilityMap = real(SusceptibilityMap).*Mask;

