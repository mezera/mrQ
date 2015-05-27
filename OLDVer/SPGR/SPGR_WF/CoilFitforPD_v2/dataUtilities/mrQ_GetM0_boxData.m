function [M01, t1, BM1,boxSize, skip , Segmask,meanVal, XX, YY, ZZ]= mrQ_GetM0_boxData(opt,T1,M0,BM,fb,smoothkernel,seg,Inclusion_Criteria)
% Load  data from the M0 and T1 file for a box cdefined by fb and opt
%
%   [M01, t1, boxSize, meanVal]= mrQ_GetM0_boxData(opt,T1,M0,fb,smoothkernel)
%
%
% Helper routine to get the data


M01=[]; t1=[]; BM1=[];boxSize=[]; meanVal=[];skip=0;Segmask=[];

%smooth (defult no )
if ~notDefined('smoothkernel') 
    if smoothkernel>0;
    for ii=1:size(M0,4)
        tmp=M0(:,:,:,ii);
        M0(:,:,:,ii)=smooth3(tmp,'g',smoothkernel);
    end
    if ~isempty(T1)
        T1(:,:,:)=smooth3(T1,'g',smoothkernel);
    end
    end
end

if notDefined('Inclusion_Criteria') 
Inclusion_Criteria=[0.8 200];
end

XX(1)=opt.X(fb(1),fb(2),fb(3))-opt.HboxS(1);
XX(2)=opt.X(fb(1),fb(2),fb(3))+opt.HboxS(1);
YY(1)=opt.Y(fb(1),fb(2),fb(3))-opt.HboxS(2);
YY(2)=opt.Y(fb(1),fb(2),fb(3))+opt.HboxS(2);
ZZ(1)=opt.Z(fb(1),fb(2),fb(3))-opt.HboxS(3);
ZZ(2)=opt.Z(fb(1),fb(2),fb(3))+opt.HboxS(3);
    
    %get the location of the box we work on in image space (x,y,z)

% Pull out the data
M01 = M0(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2),:);

% This is the 4D size  of the box
boxSize = size(M01);

% Sorting the M0 data according to the mean value of the data.
% Biggest SNR to smallest
[meanVal, coilIndex] = sort(squeeze(mean(mean(mean(M01)))),'descend');
M01     = M01(:,:,:,coilIndex);
meanVal = meanVal(coilIndex);
    if ~isempty(T1)
t1=T1(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2));
    end
BM1=logical(BM(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2)));
    

if ~notDefined('seg') 
seg=seg(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2));
Segmask=zeros(size(seg));
   TissueType= unique(seg)';
   k=1;
   for ii=TissueType
       if  length(find(seg==ii))>20 % if there is almost no voxel for a tisue type we won't use it for regularization
           Segmask(find(seg==ii))=k;
           k=k+1;
       end
   end
   if  length(find(Segmask))<(length(Segmask(:))*0.6) % we need at R1 data for regularization if we don't have it we wo't ue this box
       skip=1;
   end
   
end   
   %% if R1 values are worng we will skip the voxels
       R1=1./t1;
    Bad=isnan(R1) | isinf(R1) | R1==0;
        BM1(Bad)=0;
   


if length(find(BM1))<length(BM1(:))*Inclusion_Criteria(1) ||  length(find(BM1))<Inclusion_Criteria(2)  % not enghf voxels
    skip=1;
end
    

