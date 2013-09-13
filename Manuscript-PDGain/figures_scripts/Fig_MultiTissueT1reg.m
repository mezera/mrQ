% Fig_MultiTissueT1reg
%
% In fitting real brain data, we might have different linear relationships
% between the different types of tissue, say gray and white.  This implies
% that each would need its own separate parameters for the linear
% regularization term.
%
% Here we show how to deal with this condition and what is the cost of
% regularization with a single relationship when more then T1 PD relations
% exist.
%
% Like in case of different brain tissue types, we test what could be the
% error by simulation.
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

% This produces the key parameters for the polynomial approximations.  The
% returned variables includes the polynomial basis, pBasis, the M0 data,
% M0S_v, additional parameters, such as the box size.
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

%% Split R1 into three tissue types with different linear R1-PD relations

% In general, we could do more testing with other nonlinear relations.
%
tissue1 = find(PD<=0.7);
tissue2 = find(PD>0.7 & PD<=0.77);
tissue3 = find(PD>0.77);

R1 = zeros(size(PD));
R1(tissue1)  = (2.5./PD(tissue1)) - 2.26;
R1(tissue2)  = (1.3./PD(tissue2)) - 0.9;
R1(tissue3 ) = (1./PD(tissue3)) - 0.5;
R1=R1./1000;

mask = zeros(size(PD));
mask(tissue1)=1;
mask(tissue2)=2;
mask(tissue3)=3;

%% Simulate coil gain

% We use the poylnomial fits to the phantom data a typical coil function

% Select a set of coils
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

% Get those coil poylnomyal coeficents
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomials
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;


%% Simultae MRI SPGR signal with and with out noise
noiseLevel = 2;   % ?? Units???

% Simultate the M0 and T1 fits of multi SPGR images.
% In this case, we put R1 into the simulation so there is a built in
% nonlinear or set of distinct PD-R1 relationships.
[MR_Sim]= simSPGRs(G,PD(:),[],R1(:),[],[],noiseLevel,true);


%% Fit, under the assumption of only one tissue

% NL is a new structure with the coil coefficients (g), PD and coil image
% (G) of the volume.
NL   = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis);

%% Solve again, but add a T1 (1/R1) regularization term assuming linear

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

% mrvNewGraphWin;
% plot(PD(:)./mean(PD(:)), NL_T1reg.PD(:)./mean(NL_T1reg.PD(:)),'*')

%% Repeat the regularization accounting for the different tissue types

% This is the best we could do because here we have the proper segmentation
% of the three PD-T1 relationships.

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,[],mask,[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)))

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error');
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set
[NL_T1regTissues.PD,~,NL_T1regTissues.G,NL_T1regTissues.g, NL_T1regTissues.resnorm,NL_T1regTissues.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg),mask);


%%  reshape and scale the PD fits

% PD with T1 reg
boxSize = repmat(phantomP.rSize,1,nDims);

PD_T1reg  = reshape(NL_T1reg.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg(:));
PD_T1reg  = PD_T1reg.*scale;

PD_T1regTisuue  = reshape(NL_T1regTissues.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1regTisuue(1,1,1);
scale     = mean(PD(:)./PD_T1regTisuue(:));
PD_T1regTisuue  = PD_T1regTisuue.*scale;

PD_Noreg  = reshape(NL.PD, boxSize);
%scale     = PD(1,1,1)/PD_Noreg(1,1,1);
scale     = mean(PD(:)./PD_Noreg(:));

PD_Noreg  = PD_Noreg.*scale;
% for visualizaton i will clip it to be up to 1005 error
PD_Noreg(PD_Noreg>2)=2;
PD_Noreg(PD_Noreg<0)=0;

%% Make the figure

mrvNewGraphWin;

MM = minmax([PD_T1reg PD PD_T1regTisuue]);
hold on
plot(PD_T1regTisuue(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(PD_T1reg(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','b')

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca); xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)]);
axis image; axis square
legend('PD estimate  T1 reg','PD estimate  T1 reg Tissues','Location','NorthWest')

%%
mrvNewGraphWin;

MM = minmax([PD_T1reg PD PD_T1regTisuue]);
hold on

plot(PD_Noreg(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','r')
plot(PD_T1regTisuue(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(PD_T1reg(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','b')

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca); 
xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)]);
axis image; axis square
legend('PD estimate without T1 reg','PD estimate  T1 reg','PD estimate  T1 reg Tissues','Location','NorthWest')

%% The coefficient of determination (R^2)

% No regularization
CV1 = (calccod(PD_Noreg(:),PD(:))/100).^2

% Knowing there are three tissues
CV2 = (calccod(PD_T1regTisuue(:),PD(:))/100).^2

% Assuming only one tissue type
CV3 = (calccod(PD_T1reg(:),PD(:))/100).^2

%% END