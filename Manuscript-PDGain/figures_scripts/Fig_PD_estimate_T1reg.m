%% Figure 3
%
% Illustrate how PD can be measured from noise M0 while estimate the coil
% gain using bilinear solutions with a T1 regularization term.
%
% The problem we solve by the T1 regularization is the effect of noise on
% this solutions.
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
[PD, R1] = mrQ_simulate_PD('6',phantomP.nVoxels);



%% Simulate coil gain using the poylnomial fits to the phantom data
% These are typical coil functions


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
[PDinit, g0]=Get_PDinit(0,[],4,MR_Sim.M0SN,phantomP.pBasis);
%% Solve the bilinear  problem with no regularization

% NL is a new structure with the coil coefficients (g), PD and coil image
% (G) of the volume.
NL   = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis);


% No simulated noise
NL_noNoise = pdBiLinearFit_lsqSeach(MR_Sim.M0S,phantomP.pBasis);


%% Solve again, but add a T1 (1/R1) regularization term

%  The model we regularize asserts that there is a relationship between PD
%  and R1 (1/T1).  We believe, based on the literature, that there is a
%  relationship between PD and R1
%
%     1/PD = c1*(1/R1) + c2;
%
% The coefficients, c1 and c2, can be anything.  We regularize only on the
% idea that the linear relationship holds within each small box.  The
% coefficients can vary as we measure in different boxes.
%

% 1. Find the best weight (lambda) for the T1 reg.
%    We find this fit by X-validation
kFold   = 2; % X-validate on half the data

% Possible weights to test.  We will choose the one that cross-validates
% best. 
lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0]; 

% We are solving for c1 and c2.  If we put 1/PD in a column, P, and 1/R1 in
% column, R, and ones in a column, we have a linear equation
%
%        P = [R,Ones] * [c1,c2]'
%        P = Rmatrix * c
%

% We call [R,Ones] the R matrix, so pinv(R)*P = c
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
%%
 mrvNewGraphWin; 
  semilogy(lambda,X_valdationErr(2,:),'ko-', 'MarkerSize',10); xlabel('wieght','FontSize',16);ylabel(' Cross Validation RMSE','FontSize',16); 
axis image; axis square
xlim([-4e2 11e3])
ylim([min(X_valdationErr(2,:))*0.9 max(X_valdationErr(2,:))*1.1])
grid on
set(gca,'FontSize',16)

%%  reshape and scale the PD fits

% PD with T1 reg
PD_T1reg  = reshape(NL_T1reg.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg(:));

PD_T1reg  = PD_T1reg.*scale;

% PD with out T1 reg
PD_Noreg  = reshape(NL.PD, boxSize);
%scale     = PD(1,1,1)/PD_Noreg(1,1,1);
scale     = mean(PD(:)./PD_Noreg(:));
PD_Noreg  = PD_Noreg.*scale;

%5
PD_NoregNoNoise  = reshape(NL_noNoise.PD, boxSize);
%scale     = PD(1,1,1)/PD_Noreg(1,1,1);
scale     = mean(PD(:)./PD_NoregNoNoise(:));
PD_NoregNoNoise  = PD_NoregNoNoise.*scale;
%%  make the figure

mrvNewGraphWin;

MM = minmax([PD_T1reg PD PD_Noreg]);
hold on
plot(PD_Noreg(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')
plot(PD_T1reg(:),PD(:),'or','MarkerSize',10)

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca); xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)]);
axis image; axis square
legend('PD estimate without T1 reg','PD estimate with T1 reg','Location','NorthWest')

%%
mrvNewGraphWin
plot(PD_Noreg(:),PD(:),'ko' ,'MarkerSize',10,'MarkerFaceColor','k')
xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca);
axis image; axis square
set(gca,'FontSize',16)




%%
mrvNewGraphWin
plot(PD_T1reg(:),PD(:),'ko' ,'MarkerSize',10,'MarkerFaceColor','k')
xlabel('Estimated PD','FontSize',16); ylabel('True PD','FontSize',16);
identityLine(gca); 
axis image; axis square
set(gca,'FontSize',16)

%%
mrvNewGraphWin
plot(PD_NoregNoNoise(:),PD(:),'ko' ,'MarkerSize',10,'MarkerFaceColor','k')
xlabel('Estimated PD','FontSize',16); ylabel('True PD','FontSize',16);
identityLine(gca); 
axis image; axis square
set(gca,'FontSize',16)

%%

CV1=(calccod(PD_Noreg(:),PD(:))/100).^2

CV2=(calccod(PD_T1reg(:),PD(:))/100).^2

return

%% End

%% ridge reg (not helpfull)
% kFold=2;
% lambda1= [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0];
% [X_valdationErr,   gEstT, resnorm, FitT useX, kFold ]=pdX_valdationRidgeLoop( lambda1,kFold,MR_Sim.M0SN,phantomP.pBasis);
% %        figure;  semilogy(lambda1,X_valdationErr(2,:),'*-');         xlabel('lambda');ylabel('CV error');X_valdationErr(2,:)./min(X_valdationErr(2,:))
% 
% BestReg = find(X_valdationErr(2,:)==min(X_valdationErr(2,:)))% 
% 
% 
% NL_Ridge = pdBiLinearFit_lsqRidgeSeach(MR_Sim.M0SN,phantomP.pBasis,gEstT(:,:,1,BestReg),[],lambda1(BestReg));

%%
