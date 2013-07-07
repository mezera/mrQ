% make a toy model
%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 1;      % Second order is good for up to 5 samples
nSamples = 1;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
BasisFlag = 'qr';

printImages = false;
smoothkernel=[];
% This produces the key variables for comparing data and polynomial
% approximations. We will turn it into a function before long.
% Variables include M0S_v, pBasis, params, SZ
[OutPut] = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
% mrvNewGraphWin; imagesc(OutPut.pBasis);
% tmp = reshape(OutPut.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))
percentError = 100*OutPut.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)

%% 2) simulte M0
Par=OutPut.params(:,1:4);
%Par(1,:)=Par(1,:)./100; % what if we keep the constant close to the other values 
G=OutPut.pBasis*Par;
nVoxels=size(G,1);
nCoilsS=size(G,2);

%PD = ones(nVoxels,1);
% PD = 'single point';
% PD = 'small region';
PD = 'linear slope';
noiseLevel = 10;
[M0SN, M0S, SNR, PDsim]= simM0(G,PD,noiseLevel,true);

PDsim = reshape(PDsim,OutPut.SZ(1:3));
showMontage(PDsim);
%     M0S4D = reshape(M0S,[OutPut.SZ(1:3) nCoilsS]);
% showMontage(M0S4D(:,:,:,1));
% showMontage(M0S4D(:,:,:,2));



%% 3)fit the sulotion by bilinear solver
maxLoops = 100;
sCriterion = 1e-3;  % Stopping criterion    
Lambda = .1000;

D=diag(OutPut.W);D(1,1)=0.01;
BLSim = pdBiLinearFit(M0SN, OutPut.pBasis,...
    Lambda, maxLoops, sCriterion, [], 1, Par,D);

PDfit = reshape(BLSim.PD,OutPut.SZ(1:3));
 showMontage(PDfit);
 showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
 sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))
 RMSE=sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE)] )
