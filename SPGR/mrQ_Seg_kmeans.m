function [mrQ,seg]=mrQ_Seg_kmeans(mrQ,BMfile,CSFVal,T1file,M0file,outDir)

%%clip brain to vernrucales area then segment using k means

% 
% The function segment by FSL to three tissue take the CSF tissue restrict
% it by the T1 values the CSF is also restricted to be in the center of the
% brain in a box of ~ 60x80x40mm asuming the brain is in ACPC space this is
% where the ventrical are.
% 
% 

if notDefined('CSFVal')
    CSFVal=0.35; %R1 in -sec water in body temp  (minimum value ) any think above that will cosider as water for segmentation (the same tissue)
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

%%
seg=zeros(size(R1));

notdone=0;
mask=logical(R1); mask= mask &  R1>CSFVal & BM;
while notdone==0
    
    [IDX,C] =kmeans(R1(mask),3);
    seg(mask)=IDX;
    
    notdone=1;
    % check we don't get a strange cluster that is very small in size (if
    
    % so this might be just noise)
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

% create segmentation file of the clipped brain 
seg=zeros(size(R1));
seg(mask)=IDX;

% the tissue with the highest value is white matter, with the lowest va;ue
% is GM, and the tissue with the intermediate values is deep nuclei and
% tissue between the WM and the GM

[val,idx]=sort(C);
GMclass=1;DEEPclass=2;WMclass=3; CSFclass=4;

seg(seg==idx(1))=4; seg(seg==idx(2))=5;seg(seg==idx(3))=6;
seg(seg==4)=GMclass; seg(seg==5)=DEEPclass; seg(seg==6)=WMclass; 

%%
CSF=logical(R1); CSF= CSF&  R1<CSFVal & BM;
seg(CSF)=CSFclass; % any region that is mostly water.

segfile = fullfile(outDir,'T1w_tissue.nii.gz');
dtiWriteNiftiWrapper(single(seg), xform, segfile)
%%


%%

% clip mask size
if notDefined('boxsize')
    boxsize(1)=30;
    boxsize(2)=40;
    boxsize(3)=20;
end
sz=size(CSF); szH=round(sz./2);
XX=boxsize(1)./round(mmPerVox(1));
YY=boxsize(2)./round(mmPerVox(2));
ZZ=boxsize(3)./round(mmPerVox(3));

CSF(szH(1)+XX:end,:,:)=0;
CSF(1:szH(1)-XX,:,:)=0;

CSF(:,1:szH(2)-YY,:)=0;
CSF(:,szH(2)+YY:end,:,:)=0;

CSF(:,:,1:szH(3)-ZZ)=0;
CSF(:,:,szH(3)+ZZ:end)=0;

% smooth in space
[CSF1] = ordfilt3D(CSF,6);
CSF1=CSF &  CSF1;

CSF2= CSF1 & R1<0.25 & R1>0.2 & M0<prctile(M0(BM),99);

if notDefined('csffile')
    csffile = fullfile(outDir, 'csf_seg_T1.nii.gz');
end

dtiWriteNiftiWrapper(single(CSF2), xform, csffile);

%%

if length(find(CSF2))<200
           fprintf(['\n Warning: We could find only ' num2str(length(find(CSF1))) ' csf voxel this make the CSF WF estimation very noisy, cosider to edit csf_seg_T1.nii.gz fiel see below \n']);
end

% larger ventricals region. this mask os good for cases when CSF ventrical are hard to identify. it may reduce the accuracy 
CSF= CSF & R1<0.25 & R1>0.2 & M0<prctile(M0(BM),99);

%

csffile1 = fullfile(outDir, 'csf_seg_T1_large.nii.gz');
dtiWriteNiftiWrapper(single(CSF), xform, csffile1);


%%
mrQ.csf_large=csffile1;
mrQ.csf=csffile;
mrQ.T1w_tissue=segfile;









