%%% Sup Figure 2
%
% cheak the effectivity of ridge regration
%
% 
%
% AM/BW Vistaosft Team, 2013
%%  Make sure mrQ is on the path
addpath(genpath(fullfile(mrqRootPath)));
%% Generate example parameters for the coils from the phantom data

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box 
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations.  The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);

boxSize = repmat(phantomP.rSize,1,nDims);

%% Simulate PD 
[X,Y, Z] = meshgrid(-nSamples:nSamples,-nSamples:nSamples, -nSamples:nSamples);
R  = sqrt(X.^2 + Y.^2 + Z.^2);

% R is the distance from the center.  We make a rectified sinusoid from the
% center to the edge.  We set all the NaN values to 1.  We then take the
% sixth root to squeeze the dynamic range to be reasonable.
PD = sin(R)./R; 
PD(isnan(PD) )= 1;
PD = abs(PD);
PD = PD .^ (1/6);
% sqrt(sqrt(sqrt(PD)));

%% Simulate coil gain using the poylnomial fits to the phantom data
% These are typical coil functions

% Select a set of coils 
% Arbitrary choice: coils = [1 3 5 9]; 
% We can sort coils by minimalcorrelation between the coils to find the best set.

% We use this algorithm to
nUseCoils = 4;                         % How many coils to use
c = nchoosek(1:16,nUseCoils);          % All the potential combinations of nCoils
Cor = ones(nUseCoils,size(c,1))*100;   % Initiate the Cor to max
for kk=1:size(c,1)             % loop over coils combinations
    
    % Correlations between the the measured M0 for each of the coils in
    % this combination.
    A = (corrcoef(phantomP.M0_v(:,c(kk,:)))); 
    
    % Sum the abs correlation after correcting for number of coils and the
    % identity correlations
    Cor(nUseCoils,kk) = sum(sum(abs(triu(A) - eye(nUseCoils,nUseCoils))));
    
end

% Find the minimum correlation (min abs corr give us the set with corr that
% are closer to zero).  Choose those coils.
[v, minCor] = min(Cor(nUseCoils,:)); 
coils       = c(minCor,:);

% Get the poylnomial coeficents for those coils
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomial
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;

%% Simulate MRI SPGR signal with noise

noiseLevel = 2;   % ?? Units???

% Simulate the M0 and T1 fits of multi SPGR images. 
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,true);

% MR_Sim is a structure with multiple fields that include the simulation
% inputs MR sigunal inputs and the calculations from fitting the signal
% equation. 


%%
kFold   = 2; % X-validate on half the data

% Possible weights to test.  We will choose the one that cross-validates
% best. 
lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0]; 

% Set the fiiting loop parmeters
 maxLoops=200;
 sCriterion = 1e-3;  % Stopping criterion
%
% Fit Bilinear loop
% The bilinear search is very slow. many many iteration are needed
%
 BLFit_RidgeReg = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
                 1, maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
%
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ]=pdX_valdationLoop_RidgeReg( lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,GainPolyPar,maxLoops,sCriterion)


