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


%% Make sure mrQ is on your path

%To Set up parameters for N realistic coils
nCoils = 6;
nDims  = 3;
[params, STR] = polyGetPhantomCoef(nCoils,nDims);

%% Initialize the matrices and M0 simulation
nSamples = 10;
pOrder = 2;
[pMatrix, s] = polyCreateMatrix(nSamples,pOrder,nDims);
rSize = length(s);
nVoxels = rSize^nDims;

% clear r
% Build the simulated M0 data assuming the PD are all 1.  The parameters
% are estimated from the phantom data.
M0 = zeros(nVoxels,nCoils);
for ii=1:nCoils
    M0(:,ii)= pMatrix*params(:,ii);
end

%% Visualize the coil gains

% But I don't understand the 5 in case 3.  What is this?
figure;
switch nDims
    case 3
        for ii=1:nCoils
            im = reshape(M0(:,ii),rSize,rSize,rSize);
            subplot(nCoils,1,ii); imagesc(im(:,:,5)); axis image
        end
    case 2
        for ii=1:nCoils
            subplot(nCoils,1,ii); imagesc(reshape(M0(:,ii),rSize,rSize)); axis image
        end
    otherwise
        error('Bad number of dimensions')
end

%% Solve N coils ratio with no noise

coilList = 1:3;
est = polySolveRatio(M0(:,coilList),pMatrix);

% get the estimates
[CoilCoefErr, PDerr] = polyRatioErr(est,params(:,coilList),pMatrix);

% Without noise, the error is all zero
figure; hist(PDerr(:));

% [CoilCoefErr PDerr  estMatrix  parMatrix G M0 PD  PDspaceErr]=PolyRatioErr(est,params,pMatrix)
%%  make a noise ratio
noiseLevel = .1;
for i=1:length(coilList)
    M0Noise(:,i) = M0(:,i) + randn(size(M0,1),1)*noiseLevel;
end

% This is the signal to noise in decibels
SNR = 20*log10(mean(mean(M0(:,coilList))) /noiseLevel)

% get the gain estimates for the noisy data
est = polySolveRatio(M0Noise(:,coilList),pMatrix);

[CoilCoefErr PDerr] = polyRatioErr(est,params(:,coilList),pMatrix);
figure; hist(PDerr(:))

%
% [CoilCoefErr PDerr  estMatrix  ParMatrix G M0 PD ...
%    PDspaceErr]=PolyRatioErr(est,params,pMatrix)
%
%%

