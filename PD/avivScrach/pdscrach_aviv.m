%% Working out the multiple coil case.
%
% Adding some functions to create the polynomial matrices we need.
% See the pdScratch2 for more detail.  This one is moving a bit beyond
% that.
%
% TODO
%   1.  Realistic gains if possible
%   2.  Extend to 3D
%   3.  Add a fourth coil or even make a function for N-coils
%   4.  Summarize the error distribution maybe in PD space? Coil space?
%   5.  Virtual coil analysis
%   6.  What else?
%   7.  Make the solution with eig and stuff a function
%   8.  Make the printing out and comparison a simple function
%
%

%% First example, 2nd order, 2d
addpath('data')
addpath('avivScrach/')

%To Set up parameters for N realistic coils 
Ncoils=6;
     [params STR]=GetPantomPolyCoef(Ncoils,2)
 %[params STR]=GetPantomPolyCoef(Ncoils,dim,sample,whichCoils)
 

%% Initialize the matrices and responses
nSamples = 20;
[pMatrix, s] = polyCreateMatrix(nSamples,2,2);
rSize = length(s);


for i=1:Ncoils
r (:,i)= pMatrix*params(:,i);
end


figure;
for i=1:Ncoils
subplot(Ncoils,1,i); imagesc(reshape(r(:,i),rSize,rSize)); axis image
end
%%

% solve N coils ratio with no noise
coilslist=[1:3];
est=solveRatio(r(:,coilslist),pMatrix);
% get the estimates
 [CoilCoefErr PDerr]=PolyRatioErr(est,params(:,coilslist),pMatrix)
% [CoilCoefErr PDerr  estMatrix  ParMatrix G M0 PD  PDspaceErr]=PolyRatioErr(est,params,pMatrix)
%%  make a noise ratio 
noiseLevel = 2;
for i=1:length(coilslist)
rNoise(:,i) = r(:,i) + randn(size(r,1),1)*noiseLevel;
end

SNR=20*log10(mean(mean(r(:,coilslist))) /noiseLevel)
est=solveRatio(rNoise(:,coilslist),pMatrix);
% get the estimates

 [CoilCoefErr PDerr]=PolyRatioErr(est,params(:,coilslist),pMatrix)
% [CoilCoefErr PDerr  estMatrix  ParMatrix G M0 PD  PDspaceErr]=PolyRatioErr(est,params,pMatrix)

%%
% or should we add the nosise the the poly The ratio is noised divided!!!
% is the noise need to be scale to the signal?
for i=1:length(coilslist)
r1Noise (:,i)= pMatrix*params(:,coilslist(i)) + randn(size(r,1))*noiseLevel;
end;

