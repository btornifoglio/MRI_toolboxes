#!/usr/bin/env bash

# directory name
output_dir='D:\Analysed_data\20202205_SA_qsmtut';

# reciprocal of noise image
fslmaths $output_dir/n_std.nii.gz \
        -recip \
        $output_dir/recip_n_std.nii.gz \

# calculate mean of inverse noise map
mean_inverse_noise=$(fslstats $output_dir recip_n_std.nii.gz -M)
# mean_inverse_noise=$(fslstats $output_dir/recip_n_std.nii.gz -k $output_dir/mask.nii.gz -M)

# calculate threshold
threshold=$(awk "BEGIN {print $mean_inverse_noise/3}")

echo Mean 1/noise = $mean_inverse_noise
echo Threshold = $threshold

# calculate new threshold mask
fslmaths $output_dir recip_n_std.nii.gz -thr $threshold -bin $output_dir mask_threshold.nii.gz

# limit to slices mask_tube.nii.gz
# fslmaths $output_dir/mask_threshold.nii.gz -mul $output_dir/mask.nii.gz -bin $output_dir/mask_threshold.nii.gz

#fsleyes $output_dir *n_std.nii.gz $output_dir mask_threshold.nii.gz #$output_dir/mask.nii.gz


# red is a variable name
# blue is command
# grey is a file
# green is something ??
# orange is an operator
