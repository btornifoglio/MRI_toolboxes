% Run "qsm_analysis" in full mode to calculate QSM

% input directory
fid_dir = '\Raw_data\20201105_SA\9';

% output directory
output_dir = 'D:\Analysed_data\20202205_SA_qsmtut\qsm';

% mask directory

mask_dir = 'D:\Analysed_data\20202205_SA_qsmtut_masks\mask_threshold.nii.gz';

% run
qsm = qsm_analysis(fid_dir, 'nchannels', 1,...
                            'save_nifti', true,...
                            'RO_reorient', true,...
                            'mask_mode', false,...
                            'mask_dir', mask_dir,...
                            'output_dir',output_dir,...
                            'run_iterTik', true);
