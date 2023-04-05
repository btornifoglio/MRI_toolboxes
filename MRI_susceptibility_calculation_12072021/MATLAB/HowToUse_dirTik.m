
% Arguments
Parameters.FieldMap = FieldMap; % in ppm!
% Optional arguments
Parameters.Mask = Mask; % default = 1
Parameters.Alpha = 0.05; % delfault = 0.05
Parameters.Resolution = [1,1,1]; % default = [1,1,1]
Parameters.B0direction = [0,0,1]; % default = [0,0,1] 

Susc = dirTik(Parameters);