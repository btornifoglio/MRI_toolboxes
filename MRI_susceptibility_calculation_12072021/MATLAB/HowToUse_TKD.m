
% Arguments
Parameters.FieldMap = FieldMap; % in ppm!
% Optional arguments
Parameters.Mask = Mask; % default = 1
Parameters.Threshold = 2/3; % delfault = 2/3
Parameters.Resolution = [1,1,1]; % default = [1,1,1]
Parameters.B0direction = [0,0,1]; % default = [0,0,1] 

Susc = TKD(Parameters);