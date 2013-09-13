%%% Sup Figure 1
%
% cheak the effect of starting points on PD fit
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
  
%% intialized by the sum of squre
 [PDinit, g0]=Get_PDinit(0,[],[],MR_Sim.M0SN,phantomP.pBasis);
    %  PDinit = sqrt(sum(M0.^2,2));    % Sum of squares
    % get initial guess

% Loop over regularization weights and calculate the X-validation error

[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,g0,[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:))) 

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error'); 
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set
[NL_T1reg_In1.PD,~,NL_T1reg_In1.G,NL_T1reg_In1.g, NL_T1reg_In1.resnorm,NL_T1reg_In1.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));

% mrvNewGraphWin; 
% plot(PD(:)./mean(PD(:)), NL_T1reg.PD(:)./mean(NL_T1reg.PD(:)),'*')


%%  intilaized PD from the PD as random 

[PDinit, g0]=Get_PDinit(0,[],2,MR_Sim.M0SN,phantomP.pBasis);
 
% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,g0,[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:))) 

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error'); 
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set
[NL_T1reg_In2.PD,~,NL_T1reg_In2.G,NL_T1reg_In2.g, NL_T1reg_In2.resnorm,NL_T1reg_In2.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));
% mrvNewGraphWin; 
% plot(PD(:)./mean(PD(:)), NL_T1reg.PD(:)./mean(NL_T1reg.PD(:)),'*')

%% itiate by the T1 PD litrature relations (defult)
% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,[],[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:))) 

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error'); 
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set

[NL_T1reg_In3.PD,~,NL_T1reg_In3.G,NL_T1reg_In3.g, NL_T1reg_In3.resnorm,NL_T1reg_In3.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));

% mrvNewGraphWin; 
% plot(PD(:)./mean(PD(:)), NL_T1reg.PD(:)./mean(NL_T1reg.PD(:)),'*')



%%  reshape and scale the PD fits

% PD with T1 reg init: sum of sqr 
PD_T1reg_In1  = reshape(NL_T1reg_In1.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg_In1(:));

PD_T1reg_In1  = PD_T1reg_In1.*scale;


% PD with T1 reg  init: random
PD_T1reg_In2  = reshape(NL_T1reg_In2.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg_In2(:));

PD_T1reg_In2  = PD_T1reg_In2.*scale;

%for visualization clip to of to 100% error
PD_T1reg_In2(PD_T1reg_In2<0)=0;
PD_T1reg_In2(PD_T1reg_In2>2)=2;



% PD with T1 reg  init: T1 defult
PD_T1reg_In3  = reshape(NL_T1reg_In3.PD, boxSize);
%scale     = PD(1,1,1)/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg_In3(:));


PD_T1reg_In3  = PD_T1reg_In3.*scale;
% PD with out T1 reg

%%  make the figure

mrvNewGraphWin;

MM = minmax([PD_T1reg_In2 PD PD_T1reg_In1 PD_T1reg_In3]);
hold on
plot(PD_T1reg_In1(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')
plot(PD_T1reg_In2(:),PD(:),'or','MarkerSize',10,'MarkerFaceColor','r')
plot(PD_T1reg_In3(:),PD(:),'ok','MarkerSize',10,'MarkerFaceColor','k')
xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca); xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)]);
axis image; axis square
legend('init sum of sqr','init rand','init T1 ','Location','NorthWest')


return

%% End
