% If you have installed FSL somewhere
% other than  /usr/local/fsl/, change
% this first line accordingly.
fsldir = '/usr/local/fsl/'; 
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
setenv('FSLDIR', fsldir);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
path(path, fsldirmpath);
clear fsldir fsldirmpath;