function [MedSim, ImSim,Grad]=mrQ_QuantAnts(TargetIm,MovingIm,MovingScaleConstat)
%[MedSim, ImSim,Grad]=mrQ_QuantAnts(TargetIm,MovingIm,MovingScaleConstat)
%To check that two images agree the contrast gradient is calculated.

%example:
%TargetIm='/home/aviv.mezer/testAnts/033_SB/DifAlignTest/Align50/fitT1_GS/T1FitNLSPR_SEIR_Dat_T1.nii.gz';
%MovingIm='/home/aviv.mezer/testAnts/033_SB/DifAlignTest/Align50/Mask/WarpT1_LFit_HM.nii.gz'
%MovingScaleConstat=1000;
%maskfile='/home/aviv.mezer/testAnts/033_SB/DifAlignTest/Align50/fitT1_GS/T1FitNLSPR_SEIR_Dat_BrainMask1.nii.gz';
%
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016
%
%%
if notDefined('MovingScaleConstat'); MovingScaleConstat=1;end

Tim=readFileNifti(TargetIm);Tim=double(Tim.data);
Mim=readFileNifti(MovingIm);Mim=double(Mim.data*MovingScaleConstat);

Mask=ones(size(Tim));
Mask =Mask & Tim<3000 & Tim>100;
Mask(:,:,1)=0; Mask(:,:,end)=0; 

Tim(~Mask)=nan;
Mim(~Mask)=nan;
Im=Tim-Mim;
for i=1:size(Tim,3)
    Grad(:,:,i)=imgradient(Im(:,:,i));    
end

ImSim=Grad;
MedSim=nanmedian(abs(Grad(Mask)./Tim(Mask)));% reduce  mask


%%
% TargetIm='/home/aviv.mezer/testAnts/033_SB/Align/SEIR/fitT1_GS/T1FitNLSPR_SEIR_Dat_T1.nii.gz'
% MovingIm='/home/aviv.mezer/testAnts/033_SB/NewTests/T1AlignMask/WarpT1_LFit_HM.nii.gz'
% MovingScaleConstat=1000;
% maskfile='/home/aviv.mezer/testAnts/033_SB/Align/SEIR/fitT1_GS/mask.nii.gz'
%
%
%
% [MedSimG ImSimG,Tgrad,MgradGood]=mrQImSim(TargetIm,MovingIm,maskfile,MovingScaleConstat);
%
%
% MovingIm='/home/aviv.mezer/testAnts/033_SB/NewTests/T1Align/WarpT1_LFit_HM.nii.gz'
%
% [MedSimB ImSimB,Tgrad,MgradBad]=mrQImSim(TargetIm,MovingIm,maskfile,MovingScaleConstat);

% % RmB=impyramid(Mim(:,:,10),'reduce');
% % RmB2=impyramid(RmB,'red
