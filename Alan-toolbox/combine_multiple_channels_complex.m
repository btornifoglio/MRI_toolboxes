
function iField_combine = combine_multiple_channels_complex(iField)
% Alan Stone 31/01/2019
% This is taken from Fit_ppm_complex.m (line 70)

% combine multiple coils together, assuming the coil is the fifth dimension
iField_combine = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);
iField_combine = sqrt(abs(iField_combine)).*exp(1i*angle(iField_combine));