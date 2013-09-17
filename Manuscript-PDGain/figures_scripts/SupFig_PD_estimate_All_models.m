%% Supplementary Figure
%
% Compare PD estimation using different types of regression constraints
%
% We compare ridge, T1 and coil correlation and no regularization at all.
%
% The conclusion is that T1-regularization has the best overall properties.
% The absolute error is slightly smaller with ridge regression, but the
% spatial distribution of the error is systematic over the volume with
% ridge, and the spatial error is white noise with T1 regularization.
%
% No regularization is poor.  Coil correlation is OK, but not quite as good
% in some instances as T1 regularization.
%
% Requires vistasoft, mrQ and knkutils
%
% AM/BW Vistaosft Team, 2013

%%  Make sure mrQ is on the path
%
% addpath(genpath(fullfile(mrqRootPath)));
% addpath(genpath('/Users/wandell/Github/knkutils')); gitRemovePath;

%% Generate example parameters for the coils from the phantom data

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 3;% Which box
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations.  The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);

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

%% Simulate PD

% There are several PD spatial types.
% Type help mrQ_simulate_PD
[PD, R1] = mrQ_simulate_PD('2',phantomP.nVoxels);
 %showMontage(PD)

%% Simulate MRI SPGR signal with noise

noiseLevel = 2;   % ?? Units???

% Simulate the M0 and T1 fits of multi SPGR images.
[MR_Sim] = simSPGRs(G,PD(:),[],R1(:),[],[],noiseLevel,true);

% MR_Sim is a structure with multiple fields that include the simulation
% inputs MR sigunal inputs and the calculations from fitting the signal
% equation.

%% Initiate fit params

% The nonlinsqr fit with ridge fits
[PDinit, g0]=Get_PDinit(0,[],4,MR_Sim.M0SN,phantomP.pBasis);

options = optimset('Display','iter',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');
boxSize = repmat(phantomP.rSize,1,nDims);

%% Coil correlation regularization
coefdat = tril(corrcoef(MR_Sim.M0SN),-1);

RegWeight  = 1000;
TissueMask = logical(MR_Sim.M0SN(:,1));
XvalidationMask = logical(MR_Sim.M0SN);

[g, resnorm,dd1,exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,MR_Sim.M0SN,phantomP.pBasis,nUseCoils,RegWeight,TissueMask,coefdat,XvalidationMask),g0,[],[],options);

[PD_cor, Gn] = pdEstimate(MR_Sim.M0SN,phantomP.pBasis, g);

scale     = mean(PD(:)./PD_cor(:));
PD_cor    = PD_cor(:)*scale;
PD_cor    = reshape(PD_cor,boxSize);
CV_cor    = (calccod(PD_cor(:),PD(:))/100).^2

% For visualization clip to of to 100% error
PD_cor(PD_cor<0)=0;
PD_cor(PD_cor>2)=2;

%% Ridge regularization

maxLoops   = 200;
sCriterion = 1e-3; 
BLFit_RidgeReg = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
    1, maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));

scale      = mean(PD(:)./BLFit_RidgeReg.PD(:));
PD_ridge   = BLFit_RidgeReg.PD(:)*scale;
PD_ridge   = reshape(PD_ridge,boxSize);
CV_ridge   = (calccod(PD_ridge(:),PD(:))/100).^2

%for visualization clip to of to 100% error
PD_ridge(PD_ridge<0) = 0;
PD_ridge(PD_ridge>2) = 2;

%% R1 regularization
kFold   = 2; % X-validate on half the data

% Possible weights to test.  We will choose the one that cross-validates best.
lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0];
Rmatrix(1:phantomP.nVoxels,1) = 1;
Rmatrix(:,2) = double(MR_Sim.R1Fit);

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,g0,[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)))

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));

scale     = mean(PD(:)./NL_T1reg.PD(:));
PD_T1reg  = NL_T1reg.PD(:)*scale;
PD_T1reg  = reshape(PD_T1reg,boxSize);
CV_T1reg  = (calccod(PD_T1reg(:),PD(:))/100).^2

%for visualization clip to of to 100% error
PD_T1reg(PD_T1reg<0)=0;
PD_T1reg(PD_T1reg>2)=2;

%% no regularization
NL_NoReg  = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis,PDinit);
scale     = mean(PD(:)./NL_NoReg.PD(:));
PD_Noreg  = NL_NoReg.PD(:)*scale;
PD_Noreg  = reshape(PD_Noreg,boxSize);
CV_Noreg  = (calccod(PD_Noreg(:),PD(:))/100).^2

%for visualization clip to of to 100% error
PD_Noreg(PD_Noreg<0) = 0;
PD_Noreg(PD_Noreg>2) = 2;
%% Plot

%PD_cor; PD_ridge;  PD_T1reg; PD_Noreg;
mrvNewGraphWin
hold on
plot(PD(:),PD_Noreg(:),'.')
plot(PD(:),PD_cor(:),'.g')
plot(PD(:),PD_ridge(:),'.r')
plot(PD(:),PD_T1reg(:),'.k')

identityLine(gca);
axis image; axis square

legend( ['PD estimate without Reg R2 = ' num2str(CV_Noreg)]...
    ,['PD estimate corr Reg R2 = ' num2str(CV_cor)]...
    ,['PD estimate Ridge R2 = ' num2str(CV_ridge)],...
    ['PD estimate with T1 reg R2 = ' num2str(CV_T1reg)],...
    'Location','North')
%%

showMontage((PD_Noreg-PD)./PD);title('PD percent error  no regression')

showMontage((PD_cor-PD)./PD); title('PD percent error corr')

showMontage((PD_ridge-PD)./PD); title('PD percent error ridge regression')

showMontage((PD_T1reg-PD)./PD);title('PD percent error T1 regression')

%% Notes

%% END