#!/usr/bin/env bash

# This script is designed to be run after qsm_analsyis.m has been run in "mask_mode"

dirname='/Users/astone/tcd/data/analysis/qsm_dev.tissue_models_test/20190703_AS1/qsm_masks'

## TUBE MASK
# make tube mask
fslmaths $dirname/mag_combined.nii.gz -thr 0.15 -bin $dirname/mask-tube.nii.gz
# NOTE: manually fill top and bottom slice if needed
fslmaths $dirname/mask-tube.nii.gz -fillh $dirname/mask-tube_fillh.nii.gz
# limit tube mask to slices 65=>215
fslmaths $dirname/mask-tube.nii.gz -roi 0 -1 0 -1 64 126 0 -1 $dirname/mask-tube_short.nii.gz
# limit tube mask to slices 65=>215
fslmaths $dirname/mask-tube_fillh.nii.gz -roi 0 -1 0 -1 64 126 0 -1 $dirname/mask-tube_fillh_short.nii.gz

## BUBBLES MASK
# # remove background noise from r2s
# fslmaths $dirname/r2s.nii.gz -mul $dirname/mask-tube_short.nii.gz $dirname/r2s-tube_short.nii.gz
# # choose threshold to remove vessels and bubbles
# fslmaths $dirname/r2s-tube_short.nii.gz -thr 30 -bin $dirname/mask-bubbles.nii.gz
# # dilate, erode & fill holes
# fslmaths $dirname/mask-bubbles.nii.gz -dilM $dirname/mask-bubbles_dilM.nii.gz
# # fslmaths $dirname/mask-bubbles_dilM.nii.gz -dilM $dirname/mask-bubbles_dilM.nii.gz
# fslmaths $dirname/mask-bubbles_dilM.nii.gz -ero $dirname/mask-bubbles_ero.nii.gz
# # fslmaths $dirname/mask-bubbles_ero.nii.gz -ero $dirname/mask-bubbles_ero.nii.gz
# fslmaths $dirname/mask-bubbles_ero.nii.gz -fillh $dirname/mask-bubbles_fillh.nii.gz

## BUBBLES MASK 2
# remove background noise from r2s
fslmaths $dirname/r2s.nii.gz -mul $dirname/mask-tube_short.nii.gz $dirname/r2s-tube_short.nii.gz
# choose threshold to remove vessels and bubbles
fslmaths $dirname/r2s-tube_short.nii.gz -thr 30 -bin $dirname/mask-bubbles.nii.gz
# dilate, erode & fill holes
fslmaths $dirname/mask-bubbles.nii.gz -dilM $dirname/mask-bubbles_dilM.nii.gz
# smooth and threshold to remove some of the noise
fslmaths $dirname/mask-bubbles_dilM.nii.gz -kernel gauss 0.0497 -fmean $dirname/mask-bubbles_smooth.nii.gz
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
