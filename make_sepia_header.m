B0 = 3;                 % magnetic field strength, in Tesla
B0_dir = [0;0;1];       % main magnetic field direction, [x,y,z]
CF = 3*42.58*1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
TE = [.0046,.0092,.0138,.0184,.023,.0276]; % echo time for each GRE image, in second
delta_TE = .0046;     % echo spacing, in second
matrixSize = [224,224,28]; % image matrix size
voxelSize = [1,1,1];     % spatial resolution of the data, in mm

save sepia_header B0 B0_dir CF TE delta_TE matrixSize voxelSize