function [B1_h]=mrQ_fitL_GregB1(B1,SE_Xform,Res,s2,outDir,B1file,pixdim,s,mmPerVox,ttr,xform,degree,matlabpoolF)
% 
% [B1_h]=mrQ_fitL_GregB1(B1, SE_Xform, Res, s2, outDir, B1file, pixdim, s, ... 
%                        mmPerVox, ttr, xform, degree, matlabpoolF)
% 
%  %%% Move B1 from EPI space to SPGR. Smooth and upsample. 
% 
% 
% INPUTS:
%   (see the call in mrQ_fitT1M0)
%       B1,
%       SE_Xform,
%       Res,
%       s2,
%       outDir,
%       B1file,
%       pixdim,
%       s,
%       mmPerVox,
%       ttr,
%       xform,
%       degree,
%       matlabpoolF
% 
% OUTPUT:
%       B1_h
% 
% See Also:
%   mrQ_fitT1M0.m
% 
% 
% (C) Stanford University, VISTA Lab
% 

%#ok<*FNDSB>


%% Check INPUTS

if (~exist('degree','var') || isempty(degree)),
    degree=3;
end

if (exist('matlabpoolF','var') && matlabpoolF==1)
end


%% Fit T1 SPGR is epi space

for i=3:length(Res)
    dd(i-2).imData=Res{i}.im;
end

disp('1. Fiting lineary T1 and M0');
flipAngles = [s2(:).flipAngle];
tr = [s2(:).TR];
if(~all(tr == tr(1))), error('TR''s do not match!'); end
tr = tr(1);
[t1L,d] = relaxFitT1(cat(4,dd(:).imData),flipAngles,tr,ones(size(Res{1}.im)));

BML=logical(B1>0.3001 & B1<1.6999);
t1L(~BML)=0;


%% Fit a polynomyal guess to B1

[brainMaskSEIR, ~] = mrAnatExtractBrain(Res{1}.im, pixdim, 0.4);

tisuuemask = Res{1}.im<1700 & Res{1}.im>750 & t1L>0.4 & B1<1.6999 & B1>.30001 & brainMaskSEIR ;
tisuuemask = tisuuemask & abs(1-(Res{1}.im)./(Res{2}.im.*1000))<.5;

M = mean(B1(tisuuemask));
S = std(B1(tisuuemask));

tisuuemask    = tisuuemask & B1<(M+2*S) & B1>(M-2*S);
[tisuuemask1] = logical(ordfilt3D(tisuuemask,10));
tisuuemask1   = tisuuemask1 & tisuuemask;

[params,gains,rs] = fit3dpolynomialmodel(B1,tisuuemask1,3);
Imsz1 = size(B1);

[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),3);
B11 = reshape(Poly1*params(:),Imsz1);


%% Weights: make a confidence metric by T1 difference B1 variance and SEIR T1

% Wait=abs(1-B1./M);
M = mean(B1(tisuuemask1));
S = std(B1(tisuuemask1));

Z1 = (B1-M)./S;
Z1 = abs(Z1); Z1(Z1>4)=4;  Z1=(4-Z1)/4;% Z1(Z1<1.5)=.5;

RateT1 = Res{1}.im./t1L;
M = mean(RateT1(tisuuemask1));
S = std(RateT1(tisuuemask1));

Z1R = (RateT1-M)./S;
Z1R = abs(Z1R);  
Z1R(Z1R>4) = 4; 
Z1R = (4-Z1R)/4;

wh = 1000./Res{1}.im;
B12 = B11;
B12(tisuuemask1) = B1(tisuuemask1);

Wait = ones(Imsz1);
Wait(~tisuuemask1) = 0.3;
Wait(tisuuemask1) = Wait(tisuuemask1) + 0.3.*Z1R(tisuuemask1) +0.3.*Z1(tisuuemask1) +0.3.*wh(tisuuemask1);

[x1 y1 z1] = ind2sub(size(B1),find(tisuuemask1));

Xx(1) = max(min(x1)-10,1);
Xx(2) = min(max(x1)+10,Imsz1(1));
Yy(1) = max(min(y1)-10,1);
Yy(2) = min(max(y1)+10,Imsz1(2));
Zz(1) = max(min(z1)-1,1);
Zz(2) = min(max(z1)+1,Imsz1(3));


%% Fit local regressions

Wait(1:Xx(1),:,:)   = 0;      
Wait(Xx(2):end,:,:) = 0;

Wait(:,1:Yy(1),:)   = 0;      
Wait(:,Yy(2):end,:) = 0;

Wait(:,:,1:Zz(1))   = 0;      
Wait(:,:,Zz(2):end) = 0;

Wait(max(x1)+10:end,:,:) = 0;
Wait(:,1:min(y1)-10,:)   = 0;
Wait(:,:,1:min(z1)-10)   = 0;

[x y z]=ind2sub(size(B1),find(Wait));

% wh=tisuuemask1;

tt = ones(size(B1));
[x0 y0 z0] = ind2sub(size(B1),find(tt));

filter1 = round(20./round(pixdim));

w01 = localregression3d(x,y,z,B12(find(Wait)),(x0),(y0),(z0),[],[],filter1,Wait(find(Wait))); 
%  w0 = localregression3d(x0, y0 ,z0,B1(BM),x0,y0,z0);

tmp2 = NaN(size(B1));
tmp2(find(tt)) = w01;

tmp3 = B11;
tmp3(find(Wait)) = tmp2(find(Wait));


%% Write out the fits

B1epifile1 = fullfile(outDir,['B1_fit_lregG.nii.gz']);
dtiWriteNiftiWrapper(single(tmp2), SE_Xform, B1epifile1);

B1epifile2 = fullfile(outDir,['B1_fit_LregGx2.nii.gz']);
dtiWriteNiftiWrapper(single(tmp3), SE_Xform, B1epifile2);


%% Take B1 to the SPGR space and fill the gaps with global fits

B1_h = extractslices(s(1).imData,mmPerVox,tmp3,pixdim,ttr,1);

% Fill the area we don't have B1
[params,gains,rs] = fit3dpolynomialmodel(B1_h,(B1_h>0),4);
Imsz1 = size(B1_h);

% matlabpool close
[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),4);

B1match = reshape(Poly1*params(:),Imsz1);
B1_h(find(isnan(B1_h))) = B1match(find(isnan(B1_h)));

% B1file has to be passed in, with this commented out it will fail every
% time
% B1file=fullfile(outDir,['B1_fit_lregGx3.nii.gz']);
dtiWriteNiftiWrapper(single(B1_h), xform, B1file);


return
