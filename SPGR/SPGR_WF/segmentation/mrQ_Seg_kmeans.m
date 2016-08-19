function [mrQ]=mrQ_Seg_kmeans(mrQ,BMfile,CSFVal,T1file,M0file,outDir)
% function [mrQ,seg]=mrQ_Seg_kmeans(mrQ,BMfile,CSFVal,T1file,M0file,outDir)
%
% Clips brain to ventricles area, then segments using k-means with three
% clusters.
%
% This function uses FSL to segment into three tissues. It takes the CSF
% tissue and restricts it by the T1 values. The CSF is also restricted to
% be in the center of the brain in a box of approximately 60mm x 80mm x
% 40mm. Assuming that the brain is in AC-PC space, this is where the
% ventricles should be.
%
% ~INPUTS~
%          mrQ: The mrQ structure.
%       BMfile: The location of the Brain Mask file.
%       CSFVal: R1 for water at body temperature (in 1/sec). 
%                  Default is 0.35.
%       T1file: The location of the T1 file.
%       M0file: The location of the M0 file.
%       outDir: Directory of where the file should be saved to.
%
% ~OUTPUTS~
%          mrQ: The updated mrQ structure
%
% See also: mrQ_get_T1M0_files.m
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%


%% I. Loading and definitions
if notDefined('CSFVal')
    CSFVal=0.35; % R1 (in 1/sec) of water at body temp  (minimum value) 
    %anything above that will be considered to be water for segmentation (the same tissue)
end

if~exist('outDir')
outDir = mrQ.spgr_initDir;
end

if notDefined('BMfile')
   [~, ~ ,BMfile]=mrQ_get_T1M0_files(mrQ,0,0,1);
    if ~exist(BMfile,'file')
        BMfile = mrvSelectFile('r','Select the Brain Mask');
    end
end
BM = readFileNifti(BMfile); BM=logical(BM.data);
    
% if isfield(mrQ,'TissueFile')
% TM=readFileNifti(mrQ.TissueFile);
% TM=logical(TM.data);
% else
%     TM=BM;
% end


if notDefined('M0file')
    [~, M0file,~]=mrQ_get_T1M0_files(mrQ,0,1,0);
    if ~exist(M0file,'file')
        M0file = mrvSelectFile('r','Select the T1 file');
    end
end
M0 = readFileNifti(M0file);M0=double(M0.data);

if notDefined('T1file')
    [T1file, ~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
    if ~exist(T1file,'file')
        T1file = mrvSelectFile('r','Select the T1 file');
    end
end

disp(['Loading T1 data from ' T1file '...']);
T1 = readFileNifti(T1file);
xform = T1.qto_xyz;
mmPerVox = T1.pixdim;
T1 = double(T1.data);

R1=1./T1;

%% II. Perform k-means in a "while" loop

    fprintf('\n Performing segmentation for CSF file ...              \n');

seg=zeros(size(R1));

notdone=0;
mask=logical(R1); mask= mask &  R1>CSFVal & BM;

while notdone==0
    
    [IDX,C] =kmeans(R1(mask),3);
    seg(mask)=IDX;
    
    notdone=1;
    
    % Check we don't get a strange cluster that is very small in size.
    %(if so, this might be just noise)
    if  length(find(IDX==1))/length(IDX)<0.05
        mask=mask & seg~=1;
        notdone=0;
    end
    
    if  length(find(IDX==2))/length(IDX)<0.05
        mask=mask & seg~=2;
        notdone=0;
    end
    
    if  length(find(IDX==3))/length(IDX)<0.05
        mask=mask & seg~=3;
        notdone=0;
    end
end

% check if the clusters' means are too similar
if abs(1-C(1)/C(2))<0.1  &&  abs(1-C(1)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),1);
elseif abs(1-C(1)/C(2))<0.1
    [IDX,C] =kmeans(R1(mask),2);
elseif abs(1-C(1)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),2);
elseif abs(1-C(2)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),2);
end

%% III. Create segmentation file of the clipped brain 
seg=zeros(size(R1));
seg(mask)=IDX;

% The tissue with the highest value is white matter (WM), the tissue with
% the lowest value is gray matter (GM), and the tissue with the
% intermediate value is the deep nuclei and the tissue between the WM and
% the GM. 
%
% In some segmentations the lowest value is is air, intermediate is
% GM, and highest is WM. 
%
% Either way, the order is maintained and we get a segmentation of GM, WM
% and CSF.

[val,idx]=sort(C);
GMclass=1;DEEPclass=2;WMclass=3; CSFclass=4;

seg(seg==idx(1))=4; seg(seg==idx(2))=5;seg(seg==idx(3))=6;
seg(seg==4)=GMclass; seg(seg==5)=DEEPclass; seg(seg==6)=WMclass; 

CSF=logical(R1); CSF= CSF&  R1<CSFVal & BM;
seg(CSF)=CSFclass; % any region that is mostly water.

segfile = fullfile(outDir,'T1w_tissue.nii.gz');
dtiWriteNiftiWrapper(single(seg), xform, segfile)

%% IV save
mrQ.T1w_tissue=segfile;

%% 
% % % % % %% the CSF finding will be performed lated, in mrQ_WF
% % % % % % Clip mask size
% % % % % if notDefined('boxsize')
% % % % %     boxsize(1)=30;
% % % % %     boxsize(2)=40;
% % % % %     boxsize(3)=20;
% % % % % end
% % % % % sz=size(CSF); szH=round(sz./2);
% % % % % XX=boxsize(1)./round(mmPerVox(1));
% % % % % YY=boxsize(2)./round(mmPerVox(2));
% % % % % ZZ=boxsize(3)./round(mmPerVox(3));
% % % % % 
% % % % % CSF(szH(1)+XX:end,:,:)=0;
% % % % % CSF(1:szH(1)-XX,:,:)=0;
% % % % % 
% % % % % CSF(:,1:szH(2)-YY,:)=0;
% % % % % CSF(:,szH(2)+YY:end,:,:)=0;
% % % % % 
% % % % % CSF(:,:,1:szH(3)-ZZ)=0;
% % % % % CSF(:,:,szH(3)+ZZ:end)=0;
% % % % % 
% % % % % %% IV. Smoothe in space
% % % % % [CSF1] = ordfilt3D(CSF,6);
% % % % % CSF1=CSF &  CSF1;
% % % % % 
% % % % % CSF2= CSF1 & R1<0.25 & R1>0.2 & M0<prctile(M0(BM),99);
% % % % % 
% % % % % if notDefined('csffile')
% % % % %     csffile = fullfile(outDir, 'csf_seg_T1.nii.gz');
% % % % % end
% % % % % 
% % % % % dtiWriteNiftiWrapper(single(CSF2), xform, csffile);
% % % % % 
% % % % % %% V. Some issues
% % % % % 
% % % % % if length(find(CSF2))<200
% % % % %            fprintf(['\n Warning: We could find only ' num2str(length(find(CSF1))) ' csf voxels. This makes the CSF WF estimation very noisy. Consider editing csf_seg_T1.nii.gz file, see below \n']);
% % % % % end
% % % % % 
% % % % % % Larger ventricles region. 
% % % % % % This mask is good for cases when CSF ventricles are hard to identify. 
% % % % % % It may reduce the accuracy. 
% % % % % CSF= CSF & R1<0.25 & R1>0.2 & M0<prctile(M0(BM),99);
% % % % % 
% % % % % %
% % % % % 
% % % % % csffile1 = fullfile(outDir, 'csf_seg_T1_large.nii.gz');
% % % % % dtiWriteNiftiWrapper(single(CSF), xform, csffile1);
% % % % % 
% % % % % 
%% VI. Save
% % % mrQ.csf_large=csffile1;
% % % mrQ.csf=csffile;



