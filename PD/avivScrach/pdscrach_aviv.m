%% Working out the multiple coil case.
%
%
% TODO
%   1.  Realistic gains if possible -done X
%   2.  Extend to 3D  X
%   3.  Add a fourth coil or even make a function for N-coils  X
%   4.  Summarize the error distribution maybe in PD space? Coil space?  X
%   5.  Virtual coil analysis
%   6.  What else?
%   7.  Make the solution with eig and stuff a function  X
%   8.  Make the printing out and comparison a simple function  X
%  9.  Realistic noise
%  

%% First example, 2nd order, 2d
addpath('data')
addpath('avivScrach/')

%To Set up parameters for N realistic coils 
Ncoils=6;
     [params STR]=polyGetPantomCoef(Ncoils,3)
 %[params STR]=polyGetPantomCoef(Ncoils,dim,sample,whichCoils)
 

%% Initialize the matrices and responses
nSamples = 10;
[pMatrix, s] = polyCreateMatrix(nSamples,2,3);
rSize = length(s);

clear r
for i=1:Ncoils
M0 (:,i)= pMatrix*params(:,i);
end

%3D
figure;
for i=1:Ncoils
    im=reshape(M0(:,i),rSize,rSize,rSize);
subplot(Ncoils,1,i); imagesc(im(:,:,5)); axis image
end
%2D
figure;
for i=1:Ncoils
subplot(Ncoils,1,i); imagesc(reshape(M0(:,i),rSize,rSize)); axis image
end
%%
% solve N coils ratio with no noise
coilslist=[1:3];
est=polySolveRatio(M0(:,coilslist),pMatrix);
% get the estimates
 [CoilCoefErr PDerr]=polyRatioErr(est,params(:,coilslist),pMatrix)
% [CoilCoefErr PDerr  estMatrix  parMatrix G M0 PD  PDspaceErr]=PolyRatioErr(est,params,pMatrix)
%%  make a noise ratio 
noiseLevel = 2;
for i=1:length(coilslist)
M0Noise(:,i) = M0(:,i) + randn(size(M0,1),1)*noiseLevel;
end

SNR=20*log10(mean(mean(M0(:,coilslist))) /noiseLevel)
est=polySolveRatio(M0Noise(:,coilslist),pMatrix);
% get the estimates

 [CoilCoefErr PDerr]=polyRatioErr(est,params(:,coilslist),pMatrix)
% [CoilCoefErr PDerr  estMatrix  ParMatrix G M0 PD  PDspaceErr]=PolyRatioErr(est,params,pMatrix)

%%

