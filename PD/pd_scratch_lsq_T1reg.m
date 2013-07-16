% PD run scrip for LSQ search with T1 regolariztion

%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% 2) Run the script for the pdPolyPhantomOrder
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

%% 3) simulate M0 and R1

Par = OutPut.params(:,[1:3]);

%Par(1,:)=Par(1,:)./100; % what if we keep the constant close to the other values
G = OutPut.pBasis*Par;
nVoxels = size(G,1);
nSimulatedCoils = size(G,2);

%PD = ones(nVoxels,1);
%PD = 'single point';
%PD = 'small region';
%PD = 'linear slope';
%PD = 'tissue1';
PD = 'tissue2';  % Subset of voxels

noiseLevel = 5;   % ?? Units???

[M0SN, M0S, SNR, PDsim, mask]= simM0(G,PD,noiseLevel,true);

% This is the typical linear relationship between T1 and PD
% Could go into simM0.
R1 = (2.5./PDsim) - 0.95;

% Put this in the form of a block.
PDsim = reshape(PDsim,OutPut.SZ(1:3));

%% 4)intiate the search parameters
options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun', 1e-6,'TolX', 1e-10);
      nPolyCoef = size(OutPut.pBasis,2); 
      
%  START PD
% mean of squr
%  PDsosq = sqrt(sum(M0SN.^2,2));
%  PDinit=PDsosq;
%  PDinit=PDinit(:);

%   random
% PDinit = rand(size(PDsim(:)));
% PDinit=PDinit(:);

%   segmentaion
% PDinit=nan(size(mask));
% PDinit(find(mask==1))=1;
% PDinit=PDinit(:);

%   true solution
PDinit = PDsim(:);     
      
% get inital guess
G  = zeros(nVoxels,nSimulatedCoils);
g0 = zeros(nPolyCoef,nSimulatedCoils);

% we can be specific with what we start the rest will be zeros.
mask1 = ~isnan(PDinit);
for ii=1:nSimulatedCoils
    G(mask1,ii)  = M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) = OutPut.pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end

R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = R1(:);

%% 5) LSQ fit

% Coil list
clist = [1 2];
lambda1 = 1e4;   % Weight on T1 regularization

% Searching on the gain parameters, G.
[gEst, resnorm, dd1, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearT1reg(par, M0SN(:,clist),...
    OutPut.pBasis, nVoxels, length(clist), R1basis, lambda1),...
    double(g0(:,clist)),[],[],options);

%% 6) Visualiztion     
G = OutPut.pBasis*gEst(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0SN(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE) ' the err is : ' num2str( resnorm)] )


