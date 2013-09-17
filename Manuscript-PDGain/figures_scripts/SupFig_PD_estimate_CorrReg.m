%%% Sup Figure 3
%
% cheak the effectivity of the correlation regularization
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
 [PDinit, g0]=Get_PDinit(2,MR_Sim.R1Fit,[],MR_Sim.M0SN,phantomP.pBasis);

coefdat=tril(corrcoef(MR_Sim.M0SN),-1);
    options = optimset('Display','iter',...
        'MaxFunEvals',Inf,...
        'MaxIter',Inf,...
        'TolFun', 1e-6,...
        'TolX', 1e-10,...
        'Algorithm','levenberg-marquardt');

 
[g, resnorm,dd1,exitflag] = lsqnonlin(@(par) errlocalGainUC_v2(par,MR_Sim.M0SN,phantomP.pBasis,coefdat,nUseCoils),g0',[],[],options);


 errFitNestBiLinearCorrReg(g,M0,pBasis,nCoils,RegWeight,TissueMask,nPositions)
    [PD_cor, Gn] = pdEstimate(MR_Sim.M0SN,phantomP.pBasis, g');

scale     = mean(PD(:)./PD_cor(:));
 PD_cor       =PDn(:)*scale;
 
    CV=(calccod(PD_cor(:),PD(:))/100).^2
    
    
     
    %%
    mrvNewGraphWin
    plot(PD(:),PD_cor(:),'.')
    axis image; axis square
    identityLine(gca);
xlabel('PD SIM'),ylabel('PD fits')

%%
boxSize = repmat(phantomP.rSize,1,nDims);

    PD_cor=reshape(PD_cor,boxSize);
        showMontage((PD_cor-PD)./PD); title('PD pracent error ')
        
        
        
        %%
        
        RegWeight=1000;
        TissueMask=logical(MR_Sim.M0SN(:,1));
        XvalidationMask=logical(MR_Sim.M0SN);
        [g, resnorm,dd1,exitflag] = lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,MR_Sim.M0SN,phantomP.pBasis,nUseCoils,RegWeight,TissueMask,coefdat,XvalidationMask),g0,[],[],options);

        [PD_cor, Gn] = pdEstimate(MR_Sim.M0SN,phantomP.pBasis, g);
        
        scale     = mean(PD(:)./PD_cor(:));
        PD_cor       =PDn(:)*scale;
        
        CV=(calccod(PD_cor(:),PD(:))/100).^2