function opt=R1Seg(opt,CSFVal)
% opt=R1Seg(opt,CSFVal)
%
% Uses kmeans to segment the R1 signal of the brain.
%
% ~INPUTS~
%       opt: A structure of parameters. 
%    CSFVal: The R1 value of water at body temperature, in units of 1/sec 
%               (Default is 0.35)
%
% ~OUTPUTS~
%      opt: The updated opt structure.


if notDefined('CSFVal')
    CSFVal=0.35; %R1, in 1/sec, of water at body temp (minimum value).
                 %Anything above that value will be considered as water for 
                 %the purpose of segmentation (the same tissue)
end

BM=readFileNifti(opt.BMfile);
BM=logical(BM.data);
T1=readFileNifti(opt.T1file);

if isfield(opt,'TissueFile')
    TM=readFileNifti(opt.TissueFile);
    TM=logical(TM.data);
else
    TM=BM;
end

R1=1./T1.data;
seg=zeros(size(R1));
out=BM & ~TM;
CSF= R1<CSFVal &~out; % Any tissue with T1> ~0.285 is mostly water, and is 
                      % a different tissue (CSF) than brain tissue

R1(R1>2.5)=2.5;
mask= TM & ~CSF & R1<2 ;% In the brain, not CSF and not very low (bon or noise)

notdone=0;

%segment the R1 that is not CSF by k-means (k=3)
while notdone==0
    
    [IDX,C] =kmeans(R1(mask),3);
    seg(mask)=IDX;
    
    notdone=1;
    % Check that we don't get a strange cluster that is very small in size 
    % (If so, this might be just noise)
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

if abs(1-C(1)/C(2))<0.1  &&  abs(1-C(1)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),1);
elseif abs(1-C(1)/C(2))<0.1
    [IDX,C] =kmeans(R1(mask),2);
elseif abs(1-C(1)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),2);
elseif abs(1-C(2)/C(3))<0.1
    [IDX,C] =kmeans(R1(mask),2);
end
seg=zeros(size(R1));

seg(mask)=IDX;
% seg(CSF)=length(C)+1;
% seg(out)=length(C)+2;

filename=fullfile(opt.outDir,'R1_seg.nii.gz');
dtiWriteNiftiWrapper(single(seg),T1.qto_xyz,filename);

opt.segfile=filename;