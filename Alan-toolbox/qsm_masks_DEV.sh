#!/usr/bin/env bash
# Alan Stone (TCD) 2019
# This script is designed to be run after qsm_analsyis.m has been run in "mask_mode"
# INPUTS:
#

if [ "$#" -ne 8 ]; then
	echo '     ===> usage: qsm_masks_DEV.sh dirname thr_tube bot_slice top_slice thr_vessels_bubbles thr_vessels bubble_mask_option sigma'
	echo '     ===> e.g. qsm_masks_DEV.sh 20190703_AS1/qsm_masks 0.15 64 126 30 10 option1 0.042'
	exit
fi

## INPUTS
dirname=$1 # masks directory produced by qsm_analsyis.m has run in "mask_mode"
thr_tube=$2 # threshold chosen to extract tube and contents from nackground in magnitude combined image
bot_slice=$3 # bottom slice of mask (indexing begins at 0)
top_slice=$4 # top slice of mask (indexing begins at 0)
thr_vessels_bubbles=$5 # threshold chosen to remove vessels and bubbles in R2s image
thr_vessels=$6 # threshold chosen to localise to vessels
bubble_mask_option=$7 # option for bubble masking [option1 option2]
sigma=$8 # smoothing kernal - sigma = FWHM/2.3548. FWHM in mm

## TUBE MASK
# make tube mask
# fslmaths $dirname/mag_combined.nii.gz -thr $thr_tube -bin $dirname/mask-tube.nii.gz
# NOTE: manually fill top and bottom slice if needed
fslmaths $dirname/mask-tube.nii.gz -fillh $dirname/mask-tube_fillh.nii.gz
# limit tube mask to slices 65=>215
fslmaths $dirname/mask-tube.nii.gz -roi 0 -1 0 -1 $bot_slice $top_slice 0 -1 $dirname/mask-tube_short.nii.gz
# limit tube mask to slices 65=>215
fslmaths $dirname/mask-tube_fillh.nii.gz -roi 0 -1 0 -1 $bot_slice $top_slice 0 -1 $dirname/mask-tube_fillh_short.nii.gz


if [ $bubble_mask_option == 'option1' ]; then

    ## BUBBLES MASK - OPTION 1 ... works well for tissue model samples with bubbles
    # remove background noise from r2s
    fslmaths $dirname/r2s.nii.gz -mul $dirname/mask-tube_short.nii.gz $dirname/r2s-tube_short.nii.gz
    # choose threshold to remove vessels and bubbles
    fslmaths $dirname/r2s-tube_short.nii.gz -thr $thr_vessels_bubbles -bin $dirname/mask-bubbles.nii.gz
    # dilate, erode & fill holes
    fslmaths $dirname/mask-bubbles.nii.gz -dilM $dirname/mask-bubbles_dilM.nii.gz
    # smooth and threshold to remove some of the noise
    fslmaths $dirname/mask-bubbles_dilM.nii.gz -kernel gauss $sigma -fmean $dirname/mask-bubbles_smooth.nii.gz
    # threshold
    fslmaths $dirname/mask-bubbles_smooth.nii.gz -thr 1 $dirname/mask-bubbles_smoothT.nii.gz
    # fill holes
    fslmaths $dirname/mask-bubbles_smoothT.nii.gz -fillh $dirname/mask-bubbles+.nii.gz
    # dilate a couple of times to make sure bubble coverage is generous
    fslmaths $dirname/mask-bubbles+.nii.gz -dilM $dirname/mask-bubbles++.nii.gz
    fslmaths $dirname/mask-bubbles++.nii.gz -dilM $dirname/mask-bubbles+++.nii.gz
    # fill holes again just incase
    fslmaths $dirname/mask-bubbles+++.nii.gz -fillh $dirname/mask-bubbles++++.nii.gz
    ## KEEP VESSELS & PBS ... remove bubbles and background MASK
    fslmaths $dirname/mask-tube_short.nii.gz -sub $dirname/mask-bubbles++++.nii.gz -bin $dirname/mask-vessels+pbs.nii.gz

    ## VESSEL MASK
    # remove bubbles & backgorund noise from r2s
    fslmaths $dirname/r2s.nii.gz -mul $dirname/mask-vessels+pbs.nii.gz $dirname/r2s-vessels+pbs.nii.gz
    # threshold to localise to vessels
    fslmaths $dirname/r2s-vessels+pbs.nii.gz -thr 10 -bin $dirname/mask-vessel.nii.gz

    ## PBS / BACKGROUND MATERIAL MASK for MEDI+0
    fslmaths $dirname/mask-tube_fillh_short.nii.gz -mul $dirname/mask-vessels+pbs.nii.gz $dirname/mask-pbs.nii.gz
    fslmaths $dirname/mask-pbs.nii.gz -sub $dirname/mask-vessel.nii.gz -bin $dirname/mask-pbs.nii.gz

elif [ $bubble_mask_option == 'option2' ]; then

    ## BUBBLES MASK - OPTION 3 ... works well for plaques with bubbles and calcifications
    # choose threshold to remove calcifications and bubbles
    fslmaths $dirname/mag_combined.nii.gz -thr $thr_tube -binv $dirname/mask-bubbles-tube.nii.gz
    # constrain to tube
    fslmaths $dirname/mask-bubbles-tube.nii.gz -mul $dirname/mask-tube_fillh_short.nii.gz $dirname/mask-bubbles.nii.gz
    # dilate, erode & fill holes
    fslmaths $dirname/mask-bubbles.nii.gz -dilM $dirname/mask-bubbles_dilM.nii.gz
    # smooth and threshold to remove some of the noise
    fslmaths $dirname/mask-bubbles_dilM.nii.gz -kernel gauss $sigma -fmean $dirname/mask-bubbles_smooth.nii.gz
    # threshold
    fslmaths $dirname/mask-bubbles_smooth.nii.gz -thr 1 $dirname/mask-bubbles_smoothT.nii.gz
    # fill holes
    fslmaths $dirname/mask-bubbles_smoothT.nii.gz -fillh $dirname/mask-bubbles+.nii.gz
    # dilate a couple of times to make sure bubble coverage is generous
    fslmaths $dirname/mask-bubbles+.nii.gz -dilM $dirname/mask-bubbles++.nii.gz
    fslmaths $dirname/mask-bubbles++.nii.gz -dilM $dirname/mask-bubbles+++.nii.gz
    # fill holes again just incase
    fslmaths $dirname/mask-bubbles+++.nii.gz -fillh $dirname/mask-bubbles++++.nii.gz

    ## KEEP VESSELS & PBS ... remove bubbles and background MASK
    fslmaths $dirname/mask-tube_short.nii.gz -sub $dirname/mask-bubbles++++.nii.gz -bin $dirname/mask-vessels+pbs.nii.gz

    ## VESSEL MASK
    # tube extract r2s
    fslmaths $dirname/mask-tube_fillh.nii.gz -mul $dirname/r2s.nii.gz $dirname/r2s-tube.nii.gz
    # smooth to remove salt and pepper noise
    sigma_r2s=`echo 1 | awk "{print $sigma * 1.75}"`
    echo $sigma_r2s
    fslmaths $dirname/r2s-tube.nii.gz -kernel gauss $sigma_r2s -fmean $dirname/r2s-smooth.nii.gz
    # remove bubbles & backgorund noise from r2s
    fslmaths $dirname/r2s-smooth.nii.gz -mul $dirname/mask-vessels+pbs.nii.gz $dirname/r2s-vessels+pbs.nii.gz
    # threshold to localise to vessels
    fslmaths $dirname/r2s-vessels+pbs.nii.gz -thr $thr_vessels -bin $dirname/mask-vessel.nii.gz

    ## PBS / BACKGROUND MATERIAL MASK for MEDI+0
    fslmaths $dirname/mask-tube_fillh_short.nii.gz -mul $dirname/mask-vessels+pbs.nii.gz $dirname/mask-pbs.nii.gz
    fslmaths $dirname/mask-pbs.nii.gz -sub $dirname/mask-vessel.nii.gz -bin $dirname/mask-pbs.nii.gz

fi
