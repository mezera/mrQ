% make a toy model
%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
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
Par=OutPut.params(:,[1 2]);
%Par(1,:)=Par(1,:)./100; % what if we keep the constant close to the other values 
G=OutPut.pBasis*Par;
nVoxels=size(G,1);
nCoilsS=size(G,2);

%PD = ones(nVoxels,1);
% PD = 'single point';
% PD = 'small region';
%PD = 'linear slope';
PD = 'tissue1';
PD = 'tissue2';

noiseLevel = 5;
[M0SN, M0S, SNR, PDsim, mask]= simM0(G,PD,noiseLevel,true);

PDsim = reshape(PDsim,OutPut.SZ(1:3));
showMontage(PDsim);
%     M0S4D = reshape(M0S,[OutPut.SZ(1:3) nCoilsS]);
% showMontage(M0S4D(:,:,:,1));
% showMontage(M0S4D(:,:,:,2));

%% 3)fit the sulotion by bilinear solver
maxLoops = 100;
sCriterion = 1e-3;  % Stopping criterion    
Lambda =.0      % 1000.500000;

D = diag(OutPut.W); %D(1,1)=0.01;
% diag(D)
% Constant terme
D(1,1) = 0;
% Linear terms -
lWeight = .01;
sWeight = 1;
cWeight = 1;
D(2,2) = lWeight;  D(4,4) = lWeight; D(7,7) = lWeight;
% Quadratic terms
D(3,3) = sWeight; D(5,5)= sWeight; D(6,6) = sWeight;
% Cross products?
D(8,8) = cWeight; D(9,9) = cWeight; D(10,10) = cWeight;
HoldforCV=0; %0.4;
%PDinit = PDsim(:);
%PDinit = rand(size(PDsim(:)));
 PDinit = [];
%   coilList=(1:10);
%  Dat=M0SN(:,coilList);
%%
 PDinit=nan(size(mask));
 PDinit(find(mask==1))=1;
PDinit=PDinit(:);
 W=[];
 Lambda =.01  
%  S=sum(M0SN,2);
% W=diag(S./max(M0SN(:))); 
% W=W.^2;


 BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
    Lambda, 1, sCriterion, PDinit(:), 1, Par,D,HoldforCV,W);
PDfit = reshape(BLSim.PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE)] )

fiterr=BLSim.M0Fit(end)
maxLoops = 100;

 BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
    Lambda, maxLoops, sCriterion, PDinit(:), 1, Par,D,HoldforCV,W);
%%
close all
maxLoops=40000;
Lambda =.001      % 1000.500000;

 BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
    Lambda, maxLoops, sCriterion, BLSim.PD(:), 0, Par,D,HoldforCV,W);
BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
    Lambda, 1, sCriterion, BLSim.PD(:), 1, Par,D,HoldforCV,W);


PDfit = reshape(BLSim.PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE)] )

fiterr=BLSim.M0Fit(end)
 %%
  coilList=(1:10);
 Dat=M0SN(:,coilList);
 maxLoops = 100;

 
 
 S=sum(Dat,2);
W=diag(S./max(Dat(:))); 
W=W.^2;


 BLSim1 = pdBiLinearFit_1(Dat, OutPut.pBasis, ...
    Lambda, maxLoops, sCriterion, PDinit(:), 1, Par(:,coilList),D,HoldforCV,W);


 
 
 
 coilList=(1:7);
 Dat=M0SN(:,coilList);
 
 
 
 S=sum(Dat,2);
W=diag(S./max(Dat(:))); 
W=W.^2;

BLSim2 = pdBiLinearFit_1(Dat, OutPut.pBasis, ...
    Lambda, maxLoops, sCriterion, BLSim1.PD(:), 1, Par(:,coilList),D,HoldforCV,W);

% PDinit = BLSim.PD(:);

% To see the optimal solution
% PDinit = PDsim(:);
% BLSim = pdBiLinearFit(M0SN, OutPut.pBasis, ...
%    Lambda, maxLoops, sCriterion, PDinit(:), 1, Par,D);
%     
%[1 2 4 7]
PDfit = reshape(BLSim1.PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE)] )




%%
Lambdas=[10 5 2 1 0.5 0.1 0.01 0.001 0];
maxLoops = 1000;
sCriterion = 1e-3;  % Stopping criterion    

D = diag(OutPut.W); %D(1,1)=0.01;
% diag(D)
% Constant terme
D(1,1) = 0;
% Linear terms -
lWeight = .01;
sWeight = 1;
cWeight = 1;
D(2,2) = lWeight;  D(4,4) = lWeight; D(7,7) = lWeight;
% Quadratic terms
D(3,3) = sWeight; D(5,5)= sWeight; D(6,6) = sWeight;
% Cross products?
D(8,8) = cWeight; D(9,9) = cWeight; D(10,10) = cWeight;
HoldforCV=0;
 PDinit = [];
for ii=1:length(Lambdas)
    
Lambda =Lambdas(ii);  
BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
    Lambda, maxLoops, sCriterion, PDinit(:), 1, Par,D,HoldforCV,wightedFlag);
PDinit = BLSim.PD(:);
end








%%  is it really better then the ratio ?

% [polyRatio] = polyCreateRatio(M0SN, OutPut.pBasis);
%  estGainCoefficients = polySolveRatio(polyRatio);
% Res = polyRatioErr(estGainCoefficients, Par, OutPut.SZ(1:3), OutPut.pBasis);
% PDfit= Res.PD;
% showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
%  sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))
%  RMSE=sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
% title(['the percent error    RMSE = '   num2str(RMSE)] )


% yes much much the ratio is much worse!!!
