function [M01, t1, BM1,boxSize, skip , meanVal]= mrQ_GetM0_boxData(opt,T1,M0,BM,fb,smoothkernel)
% Load  data from the M0 and T1 file for a box cdefined by fb and opt
%
%   [M01, t1, boxSize, meanVal]= mrQ_GetM0_boxData(opt,T1,M0,fb,smoothkernel)
%
%
% Helper routine to get the data


M01=[]; t1=[]; BM1=[];boxSize=[]; meanVal=[];skip=0;

%smooth (defult no )
if ~notDefined('smoothkernel') 
    if smoothkernel>0;
    for ii=1:size(M0,4)
        tmp=M0(:,:,:,ii);
        M0(:,:,:,ii)=smooth3(tmp,'g',smoothkernel);
    end
        T1(:,:,:)=smooth3(T1,'g',smoothkernel);
    end
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

t1=T1(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2));
BM1=logical(BM(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2)));
    



if length(find(BM1))<length(BM1(:))*0.5 ||  length(find(BM1))<100  % not enghf voxels
    skip=1;
end
    

end
