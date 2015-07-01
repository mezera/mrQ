function [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir)
%       [maskepi_File] = mrQ_B1FitMask(dirAnts,AlignFile,SET1Fitfile,outDir)




%% MASK by T1 values 
% load the Res file that have the SPGR data in epi space registered by
% ants. See: mrQ_NLANTS_warp_SPGR2EPI_RB
load(AlignFile)

brainMask=Res{1}.im>0; %This are the locations were we have SEIR T1 values

% to high to low T1 & too high too low T1 disagrement btween SPGR T1
% (without B1) and SEIR.  Also area were we are not in GM WM according to
% the T1 values. 
brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))<  1.5 & Res{1}.im<2500 & Res{1}.im>350 ;
brainMask=brainMask & (Res{1}.im./(Res{2}.im*1000))>  0.5;

%% load the SEIR fiting matrix and mask according to the fit residual and outlayer

load (SET1Fitfile);
SEIRResid=ll_T1(:,:,:,4); %fit residual
a_fitparam=ll_T1(:,:,:,2); %tissue parameter form theSEIR eqation
b_fitparam=ll_T1(:,:,:,3);%tissue parameter form theSEIR eqation
clear ll_T1

%the mask need to be without outlaier in SEIR fit.
% here we remove the voxels with big residual in SEIR fit
brainMask=~isinf(SEIRResid) & SEIRResid>0 & brainMask;
brainMask_copy=brainMask;
brainMask=brainMask & SEIRResid<prctile(SEIRResid(brainMask),95) ;

%% MAASK BY registration qulaty and warp values
% We also don't have great trust to voxels that were strongly warped when
% register the spgr to epi.
% Here we load the ANTS warp maps and remove the big warping values.

% Ants may have different output ( a version difference?), we are trying to
% find the warping files

%           May be we don't need to look for the file but to have them as input


WarpAnts=fullfile(dirAnts,'WARP_SPGR_EPI_RBWarp.nii.gz');

if exist( WarpAnts,'file')
WarpAnts=readFileNifti(WarpAnts);WarpAnts=WarpAnts.data;
WarpX=squeeze(WarpAnts(:,:,:,1,1));
WarpY=squeeze(WarpAnts(:,:,:,1,2));
WarpZ=squeeze(WarpAnts(:,:,:,1,3));


else  % it seem there is a diffrent in Ants version outputs?
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

dtiWriteNiftiWrapper(single(brainMask), ResInfo.xform, maskepi_File);


%%
     