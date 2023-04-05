%% Make Magnitude and R2* images for masking

% input directory
fid_dir = '\Raw_data\20201105_SA\9';

% output directory
output_dir = 'D:\Analysed_data\20202205_SA_qsmtut\';

% run analysis in mask mode
qsm = qsm_analysis(fid_dir, 'nchannels', 1, ...
                            'save_nifti', true,...
                            'mask_mode', true,...
                            'RO_reorient', true,...
                            'output_dir', output_dir);
