function [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir,WarpAnts)
%       [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir)
%
% In this function, a brain mask is constructed in preparation for the
% construction of the B1 fit. It uses the alignfile of the SPGR data in EPI
% space, and removes undesired extreme values.
%
%     ~INPUTS~
%       dirAnts:  The location of the SPGR init file.
%     AlignFile:  The location of the SEIR-SPGR aligned file.
%   SET1Fitfile:  The location of the SEIR-EPI fit file.
%        outDir:  Where to save the file.
% 
%    ~OUTPUTS~
%  maskepi_File:  The location of the created mask.
% 
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%


%% I. MASK by T1 values 
% Load the Res file that has the SPGR data in EPI space, as registered by
% ANTs. See: mrQ_NLANTS_warp_SPGR2EPI_RB

load(AlignFile)

brainMask=Res{1}.im>0; %This contains the locations where we have SEIR T1 values

% too high too low T1 & too high too low T1 disagreement between SPGR T1
% (without B1) and SEIR.  Also area were we are not in GM or WM according to
% the T1 values. 
%% A.M and S.F check the effect of more premisive B1 mask
%OLD version
 brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))<  1.5 & Res{1}.im<2500 & Res{1}.im>350 ;
 brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))>  0.5;
% New
%brainMask=Res{2}.im>0 & Res{1}.im>0;

%BM_path=fullfile(fileparts(SET1Fitfile),'T1FitNLSPR_SEIR_Dat_BrainMask_fromR1.nii.gz');
%BM_SEIR=readFileNifti(BM_path);
%brainMask=BM_SEIR.data>0 & Res{2}.im>0 & Res{1}.im>0;
%% II. Load the SEIR fitting matrix and mask according to the fit, residuals and outliers

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

%% III. Mask by registration quality and warp values
% We also don't trust voxels that were strongly warped when we
% registered the SPGR to EPI.

% Here we load the ANTs warp maps and remove the big warping values.

% ANTs may have different output (a version difference?). We are trying to
% find the warping files.

%  Maybe we shouldn't need to look for the files, but instead have them as
%  inputs.
WarpAntsfile=[WarpAnts,'Warp.nii.gz'];
if exist(WarpAntsfile,'file')
WarpAnts=readFileNifti(WarpAntsfile);WarpAnts=WarpAnts.data;
WarpX=squeeze(WarpAnts(:,:,:,1,1));
WarpY=squeeze(WarpAnts(:,:,:,1,2));
WarpZ=squeeze(WarpAnts(:,:,:,1,3));


else  % it seems there is a difference in ANTs version outputs?
     [a,b]=fileparts(WarpAnts);[~,b]=fileparts(b);
     WarpX=fullfile(a,[b 'xvec.nii.gz']);
     if ~exist(WarpX,'file')
         b=[b,'Warp'];
     end
  WarpX=fullfile(a,[b 'xvec.nii.gz']);WarpX=readFileNifti(WarpX);WarpX=WarpX.data;
  WarpY=fullfile(a,[b 'yvec.nii.gz']);WarpY=readFileNifti(WarpY);WarpY=WarpY.data;
  WarpZ=fullfile(a,[b 'zvec.nii.gz']);WarpZ=readFileNifti(WarpZ);WarpZ=WarpZ.data;

end
    
%% IV. Mask by warp   
% Accept only between 5th and 95th percentiles for tissue parameters
% a_fitparam and b_fitparam, and accept only between 2nd and 98th
% percentiles for warp in x, y and z orientations

 brainMask=brainMask   & b_fitparam>prctile(b_fitparam(brainMask_copy),5) &a_fitparam<prctile(a_fitparam(brainMask_copy),95)& a_fitparam>prctile(a_fitparam(brainMask_copy),5) &b_fitparam<prctile(b_fitparam(brainMask_copy),95) & ...
    WarpZ< prctile(WarpZ(brainMask_copy),98) & WarpZ> prctile(WarpZ(brainMask_copy),2) & WarpY< prctile(WarpY(brainMask_copy),98) & WarpY> prctile(WarpY(brainMask_copy),2)  & WarpX< prctile(WarpX(brainMask_copy),98) & WarpX> prctile(WarpX(brainMask_copy),2)  ;


%% V. SAVE
if notDefined('outDir')
outDir=fileparts(AlignFile);
end
maskepi_File=fullfile(outDir,['maskepiF.nii.gz']);

dtiWriteNiftiWrapper(single(brainMask), Res{1}.xform, maskepi_File);


%%
     