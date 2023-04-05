function qsm = qsm_analysis(fid_dir, varargin)
% Alan Stone (TCD) 29/5/2019
% fid_dir:  directory containing BRUKER FID file
%
% Reminder:
% PATHS
% addpath(genpath('~/albin/qsm_exvivo_dev'))
% rmpath('~/albin/qsm_libraries_1')

%===================
% SETTINGS & INPUTS
%===================
default_output_dir = './qsm'; % output directory for results
default_suffix = ''; % suffix for output niftis
default_nchannels = 1; % number of channels acquired
default_run_tkd = false; % Truncated K-space Division Dipole Inversion
default_run_tkd_UCL = false; % UCL version of Truncated K-space Division Dipole Inversion algorithm (NOTE: USE THIS ONE!)
default_run_dirTik = false; % Direct Tikhonov Regularisation method (see Anita Karsa Thesis)
default_run_iterTik = false; % Iterative Tikhonov Regularisation method (see Anita Karsa Thesis)
default_run_medi = false; % Morphology Enabled Dipole Inversion
default_run_medi0 = false; % Morphology Enabled Dipole Inversion with Zero Reference using CSF (MEDI+0)
default_save_phase_mag_channels = false; % save time by skipping mag + phase nifti generation for each channel
default_save_nifti = false;
default_save_mat = false; % save a .mat file containing phase, mag and field offset with relevent acq info
default_unwrap_type = 'laplacian'; % options 'region_growing' (MEDI); 'laplacian' (MEDI)
default_bg_field_removal = 'pdf'; % options 'pdf' (MEDI); 'lbv' (MEDI)
default_freq_calc = 'complex_fit'; % 'complex_fit' or 'echo_averaging'
default_mask_dir = []; % directory containing tissue mask
default_medi0_mask_dir = []; % directory containing zero reference mask for MEDI+0
default_save_parameter_mat = false; % keep
default_mask_mode = false; % mask mode stops after magnitude & r2s nii's generated so masks can be manually made from these images
default_RO_reorient = false; % rearrange image matrix for acquisitions with coronal orientation / head-foot RO direction (Bruker)
default_oversampling_zdir = false; % account for oversampling in the z-direction. Used in 3D GRE to stop wrap-around / aliasing in z-direction
default_apply_automask = false; % apply automatically generated mask to remove noise
default_use_nechoes = []; % only use the first n echoes of ME-GRE

ip = inputParser;
addRequired(ip,'fid_dir')
addOptional(ip,'output_dir',default_output_dir);
addOptional(ip,'suffix',default_suffix);
addOptional(ip,'nchannels',default_nchannels);
addOptional(ip,'run_tkd',default_run_tkd);
addOptional(ip,'run_tkd_UCL',default_run_tkd_UCL);
addOptional(ip,'run_dirTik',default_run_dirTik);
addOptional(ip,'run_iterTik',default_run_iterTik);
addOptional(ip,'run_medi',default_run_medi);
addOptional(ip,'run_medi0',default_run_medi0);
addOptional(ip,'save_phase_mag_channels',default_save_phase_mag_channels);
addOptional(ip,'save_nifti',default_save_nifti);
addOptional(ip,'save_mat',default_save_mat);
addOptional(ip,'unwrap_type',default_unwrap_type);
addOptional(ip,'bg_field_removal',default_bg_field_removal);
addOptional(ip,'freq_calc',default_freq_calc);
addOptional(ip,'mask_dir',default_mask_dir);
addOptional(ip,'medi0_mask_dir',default_medi0_mask_dir);
addOptional(ip,'save_parameter_mat',default_save_parameter_mat);
addOptional(ip,'mask_mode',default_mask_mode);
addOptional(ip,'RO_reorient',default_RO_reorient);
addOptional(ip,'oversampling_zdir',default_oversampling_zdir);
addOptional(ip,'automask',default_apply_automask);
addOptional(ip,'use_nechoes',default_use_nechoes);

parse(ip,fid_dir,varargin{:});

inputs = ip.Results;

%=================================
% CALCULATE PHASE FROM BRUKER FID
%=================================
fprintf('\n =======================')
fprintf('\n RUNNING qsm_analysis.m')
fprintf('\n ======================= \n')

% if run in mask_mode
if inputs.mask_mode
    inputs.output_dir = sprintf('%s_masks', inputs.output_dir);
end

% check if output directory exists
if ~exist(inputs.output_dir, 'dir')
    mkdir(inputs.output_dir)
else
    fprintf('\t=> output directory exists \n')
    fprintf('\t===> looking for new dir name ')
    fprintf('... trying  ')
    while exist(inputs.output_dir, 'dir') > 1
            inputs.output_dir = sprintf('%s+', inputs.output_dir);
            fprintf('%s, ',inputs.output_dir)
    end
    fprintf('\n\t======> found one! new output dir is %s \n',inputs.output_dir)
    mkdir(inputs.output_dir)
end

fprintf('\t=> qsm_analysis set-up: \n')
disp(inputs)

% save settings to text file
fid = fopen(sprintf('%s/qsm%s.txt',inputs.output_dir,inputs.suffix), 'wt');
fprintf(fid, '%s\n', 'SETTINGS USED FOR qsm_analysis.m');
fprintf(fid, '%s\n', struct2str(inputs));
fclose(fid);

%  Load complex data and organise from Bruker fid
if (inputs.nchannels == 1)

    % single channel
    [iField voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho] = ...
        Read_Bruker_raw(inputs.fid_dir);%,inputs.oversampling_zdir); %% BT CHANGE

    if (inputs.save_parameter_mat == 0)
        % move parameter file that is output
        [status,message,messageId] = ...
            movefile('parameter.mat',sprintf('%s/parameter_ch-%d.mat', inputs.output_dir,inputs.nchannels),'f');
    else
        !rm parameter.mat
    end

elseif (inputs.nchannels > 1)

    % multiple channels
    for channelID = 1:inputs.nchannels

        fprintf('\t===> channel %d \n',channelID)

        [iField(:,:,:,:,channelID) voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho] = ...
            Read_Bruker_raw_singlecoil(inputs.fid_dir,channelID,inputs.nchannels,inputs.oversampling_zdir);

        if (inputs.save_parameter_mat)
            % move parameter file that is output
            [status,message,messageId] = ...
                movefile('parameter.mat',sprintf('%s/parameter_ch-%d.mat',inputs.output_dir,channelID),'f');
        else
            !rm parameter.mat
        end

    end

end

% reorient matrix for saggital (head-foot) Read-Out orientations
if inputs.RO_reorient
    iField = flip(fliplr(permute(iField, [2 3 1 4 5])),3);
    % reorient Affine3D and B0_dir to match reorientation of matrix
    Affine3D = [1 0 0; 0 1 0; 0 0 1];
    B0_dir = [0; 0; 1];
end

if inputs.use_nechoes > 0
    [~,~,~,total_nechoes,~] = size(iField);
    fprintf('\t=> using the first %s echoes of %s \n',num2str(inputs.use_nechoes),num2str(total_nechoes))
    iField = iField(:,:,:,[1:inputs.use_nechoes],:);
    TE = TE(1:inputs.use_nechoes);
end

% calculate magnitude and phase
iMag = abs(iField);
Phase = angle(iField);

% magnitude & phase output nifits
if inputs.save_nifti && inputs.save_phase_mag_channels

    fprintf('\t=> saving magnitude & phase nifti for each channel  \n')

    for channelID = 1:inputs.nchannels

        save_avw(iMag(:,:,:,:,channelID), sprintf('%s/mag_ch-%d%s.nii.gz',inputs.output_dir,channelID,inputs.suffix), 'd', [voxel_size delta_TE])
        save_avw(Phase(:,:,:,:,channelID), sprintf('%s/phase_ch-%d%s.nii.gz',inputs.output_dir,channelID,inputs.suffix), 'd', [voxel_size delta_TE])

    end
end

% coil & echo combined image
iMag_coil_combined = sqrt(sum(abs(iField).^2,5));
iMag_combined = sqrt(sum(abs(iMag_coil_combined).^2,4));

% magnitude & phase output nifits
if inputs.save_nifti
    save_avw(iMag_coil_combined, sprintf('%s/mag_coil_combined%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size delta_TE])
    save_avw(iMag_combined, sprintf('%s/mag_combined%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
end


%====================
% CREATE TISSUE MASK
%====================
fprintf('\t=> calculating r2* maps \n')

% combine coils ... NOTE: check that this calculation is appropriate
iField_combine = combine_multiple_channels_complex(iField);

% arlo (MEDI toolbox)
R2s = arlo(TE, abs(iField_combine));

if inputs.save_nifti
    save_avw(R2s, sprintf('%s/r2s%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
end

% check if mask exists
if exist(inputs.mask_dir, 'file') == 2
     % load existing mask
     Mask = read_avw(inputs.mask_dir);
     % save nifti
     if inputs.save_nifti
         save_avw(Mask, sprintf('%s/mask%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
     end
else
     % mask does not exist
     fprintf('\t=> warning! no mask supplied \n')
     fprintf('\t=> warning! not using a mask will cause artefacts in final susceptibility image \n')
     % make empty mask
     Mask = ones(matrix_size);
end


%===============
% PROCESS PHASE
%===============
fprintf('\t=> processing phase \n')

% calculate noise ... NOTE: this is not currently working. The iField_combine is not a true complex combination of coils
noise_level = calfieldnoise(iField_combine, Mask);
% SNR calculation NOTE: why SNR and not just signal???
if noise_level > 1 % NOTE: check put in noise_level was 0 for high res 8 channel protocol
    disp('\n===> SNR CALCULATION??? \n')
    iField = iField/noise_level; % NOTE: Check the impact of this on results
end

% calculate TE stepsize
if inputs.use_nechoes == 1
    delta_TE = TE;
else
    delta_TE = TE(2)-TE(1);
end

%================
% FREQUENCY MAPS
%================
fprintf('\t=> frequency map calculation \n')

if strcmp('echo_averaging', inputs.freq_calc)

    % Complex Fitting
    fprintf('\t===> NOT WORKING YET!!! \n')
    % fprintf('\t===> Averaging over echoes to find Freq maps \n')
    % iFreq_raw = angle(iField_combine);
    % N_std = 1./iMag.*Mask;
    %
    % iFreq = zeros(size(iField));
    % RDF = zeros(size(iField));
    %
    % for ii = 1:6
    %     [iFreq(:,:,:,ii)] = unwrapLaplacian(iFreq_raw(:,:,:,ii), matrix_size,voxel_size); %laplacian phase unwrap
    %     %         RDF(:,:,:,ii) =LBV(iFreq(:,:,:,ii),Mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3);
    %     [RDF(:,:,:,ii), ~] = PDF(iFreq(:,:,:,ii), N_std, Mask,matrix_size,voxel_size, B0_dir); % field map using pdf
    %     RDF(:,:,:,ii) = -RDF(:,:,:,ii)./(B0*g*TE(ii));
    % end
    % RDF = mean(RDF(:,:,:,1:3),4);

elseif strcmp('complex_fit', inputs.freq_calc)

    % complex fitting
    fprintf('\t===> using complex_fit for freq maps \n')

    % frequncy maps calculated using complex fitting (MEDI toolbox)
    [iFreq_raw, N_std] = Fit_ppm_complex(iField);

    if inputs.save_nifti
        save_avw(iFreq_raw, sprintf('%s/ifreq_raw%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        save_avw(N_std, sprintf('%s/n_std%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
    end

    % if run in mask mode stop here
    if inputs.mask_mode
        qsm = []; % return empty qsm image
        fprintf('\t=> qsm_analysis.m run in mask_mode \n')
        fprintf('\t=> make susceptibility mask using mag_coil_combined.nii.gz, mag_combined.nii.gz or r2s.nii.gz \n')
        fprintf('\t=> e.g fslmaths mag_combined.nii.gz -thr 4 -bin tube_mask.nii.gz \n')
        fprintf('\t=> stopping qsm_analysis.m \n')
        return
    end

    % unwrap frequency map (region growing & laplacian)
    fprintf('\t=> phase unwrapping \n')

    if strcmp(inputs.unwrap_type,'region_growing')

        fprintf('\t===> using region growing to unwrap phase (unwrapPhase.m) \n')

        % region growing (MEDI toolbox)
        iFreq = unwrapPhase(iMag_combined, iFreq_raw, matrix_size);

        if inputs.save_nifti
            save_avw(iFreq, sprintf('%s/ifreq_regiongrow%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        end

    elseif strcmp(inputs.unwrap_type,'laplacian')

        fprintf('\t===> using laplacian phase unwrap (unwrapLaplacian.m) \n')

        % laplacian phase unwrap (MEDI toolbox)
        iFreq = unwrapLaplacian(iFreq_raw, matrix_size,voxel_size);

        if inputs.save_nifti
            save_avw(iFreq, sprintf('%s/ifreq_laplacian%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        end

    end

    if inputs.automask
        fprintf('\t=> automatically masking \n')

        % Mask
        invN_std = 1./N_std;
        invN_std(invN_std == Inf) = 0;
        % NOTE: see Karsa Thesis for thresholding of masks
        % Mask 1
        aMask_1 = (invN_std) > (mean(invN_std(:)/3));
        % Mask 2
        aMask_2 = (invN_std) > mean(invN_std(:));
        % Mask 3
        aMask_3 = (invN_std) > (mean(invN_std(:)*3)); % Mask 3

        if inputs.save_nifti
            save_avw(invN_std, sprintf('%s/invN_std%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(aMask_1, sprintf('%s/AutoMask_1%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(aMask_2, sprintf('%s/AutoMask_2%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(aMask_3, sprintf('%s/AutoMask_3%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        end

    end

    % background field removal step using PDF or LBV (MEDI toolbox)
    fprintf('\t=> background field removal \n')

    if strcmp(inputs.bg_field_removal,'pdf')

        fprintf('\t===> projection onto dipole fields \n')

        % Projection onto Dipole Fields (MEDI toolbox)
        [RDF, shim, W, W_std, W_var] = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);

        if inputs.save_nifti
            save_avw(RDF, sprintf('%s/rdf_pdf%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(shim, sprintf('%s/shim_pdf%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(shim, sprintf('%s/w_pdf%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(shim, sprintf('%s/wstd_pdf%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
            save_avw(shim, sprintf('%s/w_var_pdf%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        end

    elseif strcmp(inputs.bg_field_removal,'lbv')

        fprintf('\t===> laplacian boundary value \n')

        % LBV settings
        tol = 0.01; depth = -1; peel =0; N1 = 1000; N2 = 1500; N3 = 2000;

        % LBV (MEDI toolbox)
        RDF =LBV(iFreq, Mask, matrix_size, voxel_size, tol, depth, peel, N1, N2, N3);

        if inputs.save_nifti
            save_avw(RDF, sprintf('%s/rdf_lbv%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
        end

    end

end


%=====================
% SAVE .MAT STRCUTURE
%=====================

% general info
qsm.gen_info.fid_dir = fid_dir;

% acquisition parameters
qsm.acq_info.voxel_size = voxel_size;
qsm.acq_info.matrix_size = matrix_size;
qsm.acq_info.TE = TE;
qsm.acq_info.delta_TE = delta_TE;
qsm.acq_info.Affine3D = Affine3D;
qsm.acq_info.B0_dir = B0_dir;
qsm.acq_info.TR = TR;
qsm.acq_info.NumEcho = NumEcho;

% data
qsm.data.iMag = iMag;
qsm.data.iMag_coil_combined = iMag_coil_combined;
qsm.data.iMag_combined = iMag_combined;
qsm.data.Phase = Phase;
qsm.data.iField = iField;
qsm.data.iField_combine = iField_combine;
qsm.data.R2s = R2s;
qsm.data.Mask = Mask;
qsm.data.iFreq_raw = iFreq_raw;
qsm.data.N_std = N_std;
qsm.data.iFreq = iFreq;
qsm.data.RDF = RDF;
qsm.data.shim = shim;
qsm.data.W = W;
qsm.data.W_std = W_std;
qsm.data.W_var = W_var;

if (inputs.save_mat)
    save(sprintf('%s/qsm%s.mat',inputs.output_dir,inputs.suffix),'qsm','-v7.3')
end


%==================
% DIPOLE INVERSION
%==================
fprintf('\t=> dipole inversion \n')

% check if mask exists
% ventricular CSF mask for zero referencing
% Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
% MEDI mask for zero referencing. Mask doesn't have to be CSF ...
% see: Deh et al. (2018) Magnetic Resonance in Medicine. https://doi.org/10.1002/mrm.27410

if exist(inputs.medi0_mask_dir, 'file') == 2
     % load existing mask
     Mask_CSF = read_avw(inputs.medi0_mask_dir);
     % save nifti
     if inputs.save_nifti
         save_avw(Mask_CSF, sprintf('%s/mask_medi+0%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
     end
else
     % mask does not exist
     fprintf('\t=> warning! no medi+0 mask supplied \n')
     % make empty mask
     Mask_CSF = ones(matrix_size);
end


%% TKD
if inputs.run_tkd

    fprintf('\t=> tkd \n')

    % NOTE: This step converts Fieldmap to ppm
    % This presumes that MEDI iField's & RDF in radians (see notes in OneNote Labbook)
    RDFppm = RDF./(2*pi*delta_TE*CF)*1e6; % Calculate RDF in ppm

    if inputs.save_nifti
        save_avw(RDFppm, sprintf('%s/rdf_ppm%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
    end

    if inputs.automask
        % calculation of susceptibility using TKD with Automask
        chi_automask = susceptibility_deconvolution_sim(RDFppm, Mask.*aMask_1, B0_dir, matrix_size, voxel_size);

        if inputs.save_nifti
            save_avw(chi_automask, sprintf('%s/chi_tkd_automask%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
        end
    end

    % calculation of susceptibility using TKD
    chi = susceptibility_deconvolution_sim(RDFppm, Mask, B0_dir, matrix_size, voxel_size);

    if inputs.save_nifti
        save_avw(chi, sprintf('%s/chi_tkd%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end


end


%% TKD (UCL version)
if inputs.run_tkd_UCL

    fprintf('\t=> TKD (UCL) \n')

    % NOTE: This step converts Fieldmap to ppm
    % This presumes that MEDI iField's & RDF in radians (see notes in OneNote Labbook)
    RDFppm = RDF./(2*pi*delta_TE*CF)*1e6; % Calculate RDF in ppm

    if inputs.save_nifti
        save_avw(RDFppm, sprintf('%s/rdf_ppm%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
    end

    % Arguments
    Parameters.FieldMap = RDFppm; % in ppm!
    % Optional arguments
    Parameters.Mask = Mask; % default = 1
    Parameters.Threshold = 2/3; % delfault = 2/3
    Parameters.Resolution = voxel_size; % default = [1,1,1]
    Parameters.B0direction = B0_dir'; % default = [0,0,1]

    chi_TKD = TKD(Parameters);

    if inputs.save_nifti
        save_avw(chi_TKD, sprintf('%s/chi_TKD+%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end

end


%% dirTik (Direct Tikhonov Regularised)
if inputs.run_dirTik

    fprintf('\t=> dirTik \n')

    % NOTE: This step converts Fieldmap to ppm
    % This presumes that MEDI iField's & RDF in radians (see notes in OneNote Labbook)
    RDFppm = RDF./(2*pi*delta_TE*CF)*1e6; % Calculate RDF in ppm

    if inputs.save_nifti
        save_avw(RDFppm, sprintf('%s/rdf_ppm%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
    end

    % Arguments
    Parameters.FieldMap = RDFppm; % in ppm!
    % Optional arguments
    Parameters.Mask = Mask; % default = 1
    Parameters.Alpha = 0.05; % delfault = 0.05
    Parameters.Resolution = voxel_size; % default = [1,1,1]
    Parameters.B0direction = B0_dir'; % default = [0,0,1]

    chi_dirTik = dirTik(Parameters);

    if inputs.save_nifti
        save_avw(chi_dirTik, sprintf('%s/chi_dirTik%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end

end


%% iterTik (Direct Tikhonov Regularised)
if inputs.run_iterTik

    fprintf('\t=> iterTik \n')

    % NOTE: This step converts Fieldmap to ppm
    % This presumes that MEDI iField's & RDF in radians (see notes in OneNote Labbook)
    RDFppm = RDF./(2*pi*delta_TE*CF)*1e6; % Calculate RDF in ppm

    if inputs.save_nifti
        save_avw(RDFppm, sprintf('%s/rdf_ppm%s.nii.gz',inputs.output_dir,inputs.suffix), 'd', [voxel_size])
    end

    % Arguments
    Parameters.FieldMap = RDFppm; % in ppm!
    Parameters.Mask = Mask;
    Parameters.Noise = N_std; % up to a scaling factor
    % Optional arguments
    Parameters.Alpha = 0.02; % delfault = 0.05
    Parameters.Resolution = voxel_size; % default = [1,1,1]
    Parameters.B0direction = B0_dir'; % default = [0,0,1]
    Parameters.StoppingThreshold = 0.03; % delfault = 0.03

    chi_iterTik = iterTik(Parameters);

    if inputs.save_nifti
        save_avw(chi_iterTik, sprintf('%s/chi_iterTik%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end

end

% NOTE: this is a bit of a hack
%       check above and change for consistency???
iMag = iMag_combined;

save(sprintf('%s/RDF.mat',inputs.output_dir), 'RDF', 'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'Mask', ...
     'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', '-v7.3');

%% MEDI
if inputs.run_medi

    fprintf('\t=> medi \n')

    cd(inputs.output_dir)
    % Morphology enabled dipole inversion
    [chi, cost_reg_history, cost_data_history] = MEDI_L1('lambda',1000);

    if inputs.automask
        Mask_old = Mask;
        Mask = Mask.*aMask_1;
        save(sprintf('%s/RDF.mat',inputs.output_dir), 'RDF', 'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'Mask', ...
             'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', '-v7.3');

        % Morphology enabled dipole inversion
        [chi_automask, cost_reg_history, cost_data_history] = MEDI_L1('lambda',1000);
    end


    % Optimisation MEDI
    if false
        for l1 = [500 1000 5000 10000]
            [chi_medi, cost_reg_history, cost_data_history] = MEDI_L1('lambda',l1);
            if inputs.save_nifti
                save_avw(chi_medi,sprintf('%s/chi_medi_L1_%s%s.nii.gz',...
                    inputs.output_dir,...
                    num2str(l1,'%06.f'),inputs.suffix),...
                    'd',[voxel_size])
            end
        end
    end

    cd ..

    if inputs.save_nifti
        save_avw(chi, sprintf('%s/chi_medi%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
        if inputs.automask
            save_avw(chi_automask, sprintf('%s/chi_medi_automask%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
        end
        % save_avw(chi100, sprintf('%s/CHI_MEDI%s-L100.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
        % save_avw(chi400, sprintf('%s/CHI_MEDI%s-L400.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
        % save_avw(chi1000, sprintf('%s/CHI_MEDI%s-L1000.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end
end


%% MEDI+0
if inputs.run_medi0

    fprintf('\t=> medi+0 \n')

    % ventricular CSF mask for zero referencing
    % Mask_CSF = extract_CSF(R2s, Mask, voxel_size);

    % MEDI mask for zero referencing. Mask doesn't have to be CSF ...
    % see: Deh et al. (2018) Magnetic Resonance in Medicine. https://doi.org/10.1002/mrm.27410
    Mask_CSF = read_avw(inputs.medi0_mask_dir);

    if inputs.save_nifti
        save_avw(Mask_CSF, sprintf('%s/mask_medi+0%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end

    % NOTE: this is a bit of a hack
    %       check above and change for consistency???
    iMag = iMag_combined;

    save(sprintf('%s/RDF.mat',inputs.output_dir), 'RDF', 'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'Mask', ...
         'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF', '-v7.3');

    cd(inputs.output_dir)

    % Morphology enabled dipole inversion with zero reference using CSF (MEDI+0)
    [chi, cost_reg_history, cost_data_history] = MEDI_L1('lambda',1000,'lambda_CSF',100);

    if inputs.save_nifti
        save_avw(chi, sprintf('%s/chi_medi+0%s.nii.gz', inputs.output_dir, inputs.suffix), 'd', [voxel_size])
    end

    % Optimisation MEDI+0
    if false
        for l1 = 10000
            for l2 = [1000 10000 100000]
                [chi_medi0, cost_reg_history, cost_data_history] = MEDI_L1('lambda',l1,'lambda_CSF',l2);
                if inputs.save_nifti
                    save_avw(chi_medi0,sprintf('%s/chi_medi+0_L1_%s_L2_%s%s.nii.gz',...
                        inputs.output_dir,...
                        num2str(l1,'%06.f'),num2str(l2,'%06.f'),inputs.suffix),...
                        'd',[voxel_size])
                end
            end
        end
    end

    cd ..
end

% save(sprintf('%s/QSM%s.mat', inputs.output_dir, inputs.suffix),'-v7.3')
%}