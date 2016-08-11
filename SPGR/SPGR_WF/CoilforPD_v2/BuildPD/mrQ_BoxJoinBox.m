function [PD_fit, opt]= mrQ_BoxJoinBox(Boxes,Cbox,opt,BMfile)
% [PD_fit, opt]= mrQ_BoxJoinBox(Boxes,Cbox,opt,BMfile)
%
% This is Step 4 of 6 (including Step 0) in the pipeline to build the WF
% (water fraction) map. In this step, the boxes from the previous steps
% (mrQ_boxScaleGlobLinear) are joined together to form the PD image.
%
%
% ~INPUTS~
%           Boxes:
%            Cbox:
%             opt:   The "opt" structure of optimized parameters
%          BMfile:   The location of the Brain Mask file
%
% ~OUTPUTS~
%             opt:   The updated opt structure of optimized parameters
%               G:   Location of the newly constructed Gain image
%              PD:   Location of the newly constructed PD image
%
%
% See also: mrQ_buildPD_ver2
%           Step_0: none
%           Step_1: mrQ_CalBoxPD_step1a
%           Step_2: mrQ_ScaleBoxes_step2
%           Step_4: mrQ_smoothGain_step4b
%           Step_5: mrQ_PD2WF_step5
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015


%%  

% Find the boxes that we have estimations for
wh=find(Cbox);

if notDefined('BMfile')
    BMfile=opt.BMfile;
end

  BM=readFileNifti(BMfile);  
  xform=BM.qto_xyz;

BM=BM.data;

SZ=size(BM);

% The boxes that we will work on.
% M0 matrix (voxels, max number of PD estimation of each voxel)
M0=zeros(length(BM(:)),32);
clear BM;

% Loop over the boxes and add the value for each voxel. 
% We have kept the multiple voxel values

for ii=wh'
    % location of each voxel N copy
    loc= sum(logical(M0(Boxes(ii).loc,:)),2)+1;
    
    % ind are the locations of the voxels Boxes(ii).loc(:) 
    % and the voxel copy in the M0 2D matrix
    ind= sub2ind(size(M0),Boxes(ii).loc(:),loc(:));
    
    % Add the scale box PD values to the locations in the matrix M0
    M0(ind)=Boxes(ii).PD*Cbox(ii);   
end

%  Fill the empty spots (zeors) with NaN. 
% (Not all voxels end up having the same number of estimations)
M0(find(M0==0))=nan;
%  Take the PD median value for each voxel
PD_fit=nanmedian(M0,2);

% Reshape PD from a vector of voxels into a 3D image.
PD_fit=reshape(PD_fit,SZ);

% Save
PD_filename=fullfile(opt.outDir, 'PD_Partical.nii.gz');
opt.PD_Partical=PD_filename;

dtiWriteNiftiWrapper(single(PD_fit),xform,PD_filename);

save(opt.logname,'opt')
