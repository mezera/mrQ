function [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir)
%       [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir)
%
%      INPUTS
%       dirAnts:  The location of the SPGR init file.
%     AlignFile:  The location of the SEIR-SPGR aligned file.
%   SET1Fitfile:  The location of the SEIR-EPI fit file.
%        outDir:  Where to save the file.
% 
%      OUTPUTS
%  maskepi_File:  The location of the created mask.
% 
%

%% MASK by T1 values 
% Load the Res file that has the SPGR data in EPI space, as registered by
% ANTS. See: mrQ_NLANTS_warp_SPGR2EPI_RB

load(AlignFile)

brainMask=Res{1}.im>0; %This contains the locations where we have SEIR T1 values

% too high too low T1 & too high too low T1 disagreement between SPGR T1
% (without B1) and SEIR.  Also area were we are not in GM or WM according to
% the T1 values. 
brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))<  1.5 & Res{1}.im<2500 & Res{1}.im>350 ;
brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))>  0.5;

%% load the SEIR fitting matrix and mask according to the fit, residuals and outliers

load (SET1Fitfile);
SEIRResid=ll_T1(:,:,:,4); %fit residual
a_fitparam=ll_T1(:,:,:,2); %tissue parameter from the SEIR equation
b_fitparam=ll_T1(:,:,:,3);%tissue parameter from the SEIR equation
clear ll_T1

% The mask needs to be without outliers in SEIR fit.

% Here we remove the voxels with big residuals in the SEIR fit
brainMask=~isinf(SEIRResid) & SEIRResid>0 & brainMask;
brainMask_copy=brainMask;
brainMask=brainMask & SEIRResid<prctile(SEIRResid(brainMask),95) ;

%% Mask by registration quality and warp values
% We also don't trust voxels that were strongly warped when we
% registered the SPGR to EPI.

% Here we load the ANTS warp maps and remove the big warping values.

% ANTS may have different output (a version difference?). We are trying to
% find the warping files.

%  Maybe we shouldn't need to look for the files, but instead have them as
%  inputs.

WarpAnts=fullfile(dirAnts,'WARP_SPGR_EPI_RBWarp.nii.gz');

if exist( WarpAnts,'file')
WarpAnts=readFileNifti(WarpAnts);WarpAnts=WarpAnts.data;
WarpX=squeeze(WarpAnts(:,:,:,1,1));
WarpY=squeeze(WarpAnts(:,:,:,1,2));
WarpZ=squeeze(WarpAnts(:,:,:,1,3));


else  % it seems there is a difference in ANTS version outputs?
  WarpX=fullfile(dirAnts,'WARP_SPGR_EPI_RBWarpxvec.nii.gz');WarpX=readFileNifti(WarpX);WarpX=WarpX.data;
  WarpY=fullfile(dirAnts,'WARP_SPGR_EPI_RBWarpyvec.nii.gz');WarpY=readFileNifti(WarpY);WarpY=WarpY.data;
  WarpZ=fullfile(dirAnts,'WARP_SPGR_EPI_RBWarpzvec.nii.gz');WarpZ=readFileNifti(WarpZ);WarpZ=WarpZ.data;

end
    
%% mask by warp   
brainMask=brainMask   & b_fitparam>prctile(b_fitparam(brainMask_copy),5) &a_fitparam<prctile(a_fitparam(brainMask_copy),95)& a_fitparam>prctile(a_fitparam(brainMask_copy),5) &b_fitparam<prctile(b_fitparam(brainMask_copy),95) & ...
    WarpZ< prctile(WarpZ(brainMask_copy),98) & WarpZ> prctile(WarpZ(brainMask_copy),2) & WarpY< prctile(WarpY(brainMask_copy),98) & WarpY> prctile(WarpY(brainMask_copy),2)  & WarpX< prctile(WarpX(brainMask_copy),98) & WarpX> prctile(WarpX(brainMask_copy),2)  ;

%% SAVE
if notDefined('outDir')
outDir=fileparts(AlignFile);
end
maskepi_File=fullfile(outDir,['maskepiF.nii.gz']);

dtiWriteNiftiWrapper(single(brainMask), Res{1}.xform, maskepi_File);


%%
     