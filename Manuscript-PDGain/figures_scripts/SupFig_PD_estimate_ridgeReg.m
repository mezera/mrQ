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


%% Simulate PD

% This function simulates different types of PD and R1 volumes.
% 6 is a good selection
[PD, R1] = mrQ_simulate_PD('6',phantomP.nVoxels);

%% Simulate coil gain using the poylnomial fits to the phantom data


nUseCoils = 4;                         % How many coils to use
MaxcoilNum=16;                     %last coil to consider
% These are typical coil functions
% We can sort coils by minimalcorrelation between the coils to find the best set.
% We use this algorithm to select the coils
% Find the minimum correlation (min abs corr give us the set with corr that
% are closer to zero).  Choose those coils.

coils=mrQ_select_coilsMinCorrelation(nUseCoils,MaxcoilNum,phantomP.M0_v);



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

%% Initiate fit params

% The nonlinsqr fit with ridge fits
[PDinit, g0]=Get_PDinit(0,[],1,MR_Sim.M0SN,phantomP.pBasis);

options = optimset('Display','off',...
    'MaxFunEvals',Inf,...
    'MaxIter',200,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');
boxSize = repmat(phantomP.rSize,1,nDims);

%%

maxLoops   = 200;
sCriterion = 1e-3; 
BLFit_RidgeReg = pdBiLinearFit(MR_Sim.M0SN, phantomP.pBasis, ...
    1, maxLoops, sCriterion, [], 0 ,GainPolyPar,PD(:));

scale      = mean(PD(:)./BLFit_RidgeReg.PD(:));
PD_ridge   = BLFit_RidgeReg.PD(:)*scale;
PD_ridge  = reshape(PD_ridge,boxSize);

figure;plot(PD(:),PD_ridge(:),'*');identityLine(gca);axis image; axis square

%%
fit_method_lsq=1;
BLFit_RidgeFit = pdBiLinearFit_lsqRidgeSeach(MR_Sim.M0SN,phantomP.pBasis,g0,[],1e14,fit_method_lsq,options);
scale      = mean(PD(:)./BLFit_RidgeFit.PD(:));
PD_ridgeS   = BLFit_RidgeFit.PD(:)*scale;
PD_ridgeS  = reshape(PD_ridgeS,boxSize);
figure;plot(PD(:),PD_ridgeS(:),'*');identityLine(gca);axis image; axis square

%%
kFold   = 2; % X-validate on half the data

% Possible weights to test.  We will choose the one that cross-validates
% best.
lambda = [1e20 1e16 1e14 1e10 1e8  1e4  1e0  0];

% Set the fiiting loop parmeters

[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ]=pdX_valdationLoop_Ridge_lsq( lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,g0,[],options)

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
%  BLFit_RidgeReg = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
%                  1, maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
%
[X_valdationErr,  BLFit_RidgeReg , FitT, useX, kFold ]=pdX_valdationLoop_RidgeReg( lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,GainPolyPar,maxLoops,sCriterion);
mrvNewGraphWin;semilogy(lambda,X_valdationErr(1,:),'*-'); xlabel('lambda');ylabel('X-V error');

% reidge regretion is not X validate
%%




%%
Rmatrix(1:phantomP.nVoxels,1) = 1;
% Sometimes it is single, when from NIFTI.
Rmatrix(:,2) = double(MR_Sim.R1Fit);

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,[],[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)))

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error');
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));




%% yet ridge regression is much more accurate then without it
%
%with regularization
BLFit_RidgeReg = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
    1, maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
%with out
maxLoops=1000;
BLFit_NoRidgeReg = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
    0, maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));


%%  Scale by the mean because that is always uncertain

%ridge
scale     = mean(PD(:)./BLFit_RidgeReg.PD(:));
PD_ridgre       =BLFit_RidgeReg.PD(:)*scale;

scale     = mean(PD(:)./BLFit_NoRidgeReg.PD(:));
PD_noReg       =BLFit_NoRidgeReg.PD(:)*scale;

scale     = mean(PD(:)./NL_T1reg.PD(:));
PD_T1ref       =NL_T1reg.PD(:)*scale;
%% Show the data

mrvNewGraphWin
plot(PD(:),PD_ridgre(:),'ob')


hold on
plot(PD(:),PD_T1ref(:),'ok')

% plot(PD(:),PD_noReg(:),'or')

identityLine(gca);
axis image; axis square
legend('PD estimate Ridge','PD estimate with T1 reg','PD estimate without Reg','Location','North')

%% Ridge reg has the best accuracy, better than T1 reg 

% ridge reg
% CV=(calccod(PD_ridgre(:),PD(:))/100).^2
% % no reg
% CV=(calccod(PD_noReg(:),PD(:))/100).^2
% % T1 reg
% CV=(calccod(PD_T1ref(:),PD(:))/100).^2
% T1 reg is also very good. T1 the regulaizer is also noise. that might be why it have a bit more noise then T1 reg.
% note that the error is bigger when the M0 PD fit is with high unsertenty
% (low signal) like in the CSF.

%% The ridge regression solution has some bias in space
% T1 reg is free of that.

boxSize = repmat(phantomP.rSize,1,nDims);

PD_ridgre=reshape(PD_ridgre,boxSize);
PD_noReg=reshape(PD_noReg,boxSize);
PD_T1ref=reshape(PD_T1ref,boxSize);

% with  regularization end up with error in space (in the shape of
% PD)
showMontage((PD_ridgre - PD)./PD); title('PD percent error ridge regression')

% no regularization end up with error in space (like a coil function)
showMontage((PD_noReg-PD)./PD);title('PD percent error bilinear no regression')

% no regularization end up with error in space (like a coil function)
showMontage((PD_T1ref-PD)./PD);title('PD percent error T1 regression')

%% Notes



%% END