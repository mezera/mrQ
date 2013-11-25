function [PD_fit]= mrQ_BoxJoinBox(Boxes,Cbox,opt)
%[PD_fit]= mrQ_BoxJoinBox(Boxes,Cbox,SHub,opt);
% % join the boxes (Boxes) acording to the box Costats (Cbox) calculate by
% mrQ_boxScaleGlobLinear.m .
% make a PD image from all those BOXES PD PD_FIT.

%the box that we have estimations for
wh=find(Cbox);

%brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;
SZ=size(BM);

% the box we will work on
% M0 matrix (voxels,max number of PD estimation of each voxel)
M0=zeros(length(BM(:)),32);
clear BM;
%loop over the box and add the value for each voxel. we have keep the multipal voxels values
for ii=wh'

    % loc of each voxel N copy
loc= sum(logical(M0(Boxes(ii).loc,:)),2)+1;

%ind are the locations of the voxels Boxes(ii).loc(:) and the voxel copy in
%the M0 2D matrix
ind= sub2ind(size(M0),Boxes(ii).loc(:),loc(:));

% add the scale box PD values to the locations in matrix M0  
   M0(ind)=Boxes(ii).PD*Cbox(ii);
   
end

% feel the empty spots (zeors) with nan. (not all voxel end having the same number
% of estimations ).
M0(find(M0==0))=nan;
% take the PD median value for each voxel
PD_fit=nanmedian(M0,2);

%reshape PD from a vector of voxels to a 3D image.
PD_fit=reshape(PD_fit,SZ);

%Done


