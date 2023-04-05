
% Arguments
Parameters.FieldMap = FieldMap; % in ppm!
Parameters.Mask = Mask;
Parameters.Noise = Noise; % up to a scaling factor

% The iterative fitting is weighted by the inverse of Noise to assign lower 
% weights to less reliable (noisy) voxels. 
% For multi-echo data:
% Similarly to the MEDI method, this could be e.g. the dp1 output of the 
% Fit_ppm_complex (or Fit_ppm_complex_TE) function from the MEDI toolbox. 
% dp1 is actually equivalent to eq. 14 from Kressler et al. IEEE TMI 29.2
% (2009) and can be calculated using the multi-echo magnitude images.
% For single-echo data:
% The inverse of the magnitude image can be used as the noise map input as
% voxels with noisy phase also tend to have low magnitude.

% Optional arguments
Parameters.Alpha = 0.05; % delfault = 0.05
Parameters.Resolution = [1,1,1]; % default = [1,1,1]
Parameters.B0direction = [0,0,1]; % default = [0,0,1]
Parameters.StoppingThreshold = 0.03; % delfault = 0.03

Susc = iterTik(Parameters);