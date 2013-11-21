function opt=R1Seg(opt,CSFVal)



if notDefined('CSFVal')
    CSFVal=0.35; %R1 in -sec water in boday temp  (minimum value ) any think above that will cosider as water for segmentation (the same tissue)
end

BM=readFileNifti(opt.BMfile);
BM=logical(BM.data);
T1=readFileNifti(opt.T1file);


R1=1./T1.data;
seg=zeros(size(R1));
CSF= R1<CSFVal; %any tisuue with T1> ~0.2.85 is mostlly water and is a differnt tissue (CSF) then brain tissue

R1(R1>2.5)=2.5;
mask= BM & ~CSF & R1<2 ;% in the brain not CSF and not very low (bon or noise)
notdone=0;

%segment the R1 that is not CSF by kmeans (k=3)
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
seg(CSF)=length(C)+1;
filename=fullfile(opt.outDir,'R1_seg.nii.gz');
dtiWriteNiftiWrapper(single(seg),T1.qto_xyz,filename);

opt.segfile=filename;