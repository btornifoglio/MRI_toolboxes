function SusceptibilityMap = iterTik(Parameters)
%DESCRIPTION: SusceptibilityMap=ak_Tikhonov_iter(FrequencyMap, Noise, Parameters)
%             Calculates the susceptibility map by solving a least squares
%             problem with Tikhonov regularisation using the conjugate gradient method:
%             x = argmin_x || WM(B - B0(x*d)) ||_2^2 + alpha||Mx||_2^2 where * denotes
%             convolution
%
%INPUTS:
%   Parameters(struct): Parameters.FieldMap(double array) - input image in ppm
%                       Parameters.Mask(array) - binary tissue mask
%                       Parameters.Noise(array) - noise map
%                       Parameters.Alpha(double) - regularisation parameter (default = 0.05)
%                       Parameters.Resolution(double vector) - image resolution vector (dx,dy,dz) in mm (default = 1 mm isotropic)
%                       Parameters.B0direction(double vector) - 3-element unit vector aligned with B0 (default = [0,0,1]) 
%                       Parameters.StoppingThreshold(double) - stopping threshold (norm(res)/norm(b)) for the CG method (default = 0.03) 
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
    warndlg('Please specify Parameters.Mask!','Warning')
    return;
end
if isfield(Parameters,'Noise')
    Noise = double(Parameters.Noise);
else
    warndlg('Please specify Parameters.Noise!','Warning')
    return;
end
if isfield(Parameters,'Alpha')
    Alpha = Parameters.Alpha;
else
    Alpha = 0.05;
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
if isfield(Parameters,'StoppingThreshold')
    StoppingThreshold = Parameters.StoppingThreshold;
else
    StoppingThreshold = 0.03;
end

% Calculate fitting Weights
Noise = Noise.*Mask;
Weights = 1./Noise;
Weights(isnan(Weights)) = 0;
Weights(isinf(Weights)) = 0;
Weights = Weights/mean(Weights(Mask==1));

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

% Calculate the right-hand side of the equations: b = DW^2 f
b = real(fftshift(ifftn(ifftshift(kSpaceDipole.*fftshift(fftn(ifftshift(Weights.^2.*FieldMap)))))));
b = b(:);

% Define system matrix: A = DW^2D + alpha*M^2 and solve using the conjugate gradient method
A=@(xx)(SystemMatrix(kSpaceDipole,Mask,Weights,Alpha,xx));
SusceptibilityMap = ConjugateGradient(A,b,StoppingThreshold);
SusceptibilityMap = reshape(SusceptibilityMap,size(FieldMap));

% PSF correction
PSF = ifftn(ifftshift(kSpaceDipole.^2./(kSpaceDipole.^2 + Alpha)));
SusceptibilityMap = SusceptibilityMap/PSF(1).*Mask;

end

function b = SystemMatrix(kSpaceDipole,Mask,Weights,Alpha,x)

x = reshape(x,size(kSpaceDipole));
b = real(fftshift(ifftn(ifftshift(kSpaceDipole.*fftshift(fftn(ifftshift(x)))))));
b = real(fftshift(ifftn(ifftshift(kSpaceDipole.*fftshift(fftn(ifftshift(Weights.^2.*b)))))));
b = b + Alpha*Mask.*x;
b = b(:);

end

function x = ConjugateGradient(A,b,tol)

r = b;
p = b;
x = zeros(size(r));

while sqrt((r'*r)/(b'*b)) > tol
    
    Ap = A(p);
    pAp = p'*Ap;
    alpha = p'*r/pAp;
    
    x = x + alpha*p;
    r = r - alpha*Ap;
    
    
    Ar = A(r);
    beta = p'*Ar/pAp;
    p = r - beta*p;
    
end

end