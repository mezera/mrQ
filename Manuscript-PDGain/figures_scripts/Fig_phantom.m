%%% Figure 4
%
% Illustrate the poyinomyal order of coils with phantom data
%
% AM/BW Vistaosft Team, 2013



%%  Make sure mrQ is on the path
addpath(genpath(fullfile(mrqRootPath)));
%% Generate example parameters for the coils from the phantom data

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 3;% Which box 
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.


%% get the polynomyails error


for nSamples=2:10 % The box is -nSamples:nSamples
for pOrder=1:3 %  polynomial order

phantomP(nSamples,pOrder) = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
end
end

%% get the polynomyails error

for nSamples=2:10 % The box is -nSamples:nSamples
for pOrder=1:3 %  polynomial order
    
    Voulume(nSamples,pOrder)=((nSamples*2+1)*2)^3; % the box we used are voulume with -nSamples:nSamples = (nSamples*2+1) voxel  a side.
    % the resultion of the phantom scan is 2mm. so multipal by 2. and ^3
    % for voulume
PE(nSamples,pOrder)=phantomP(nSamples,pOrder).percentError;
end
end


mrvNewGraphWin;
hold on


plot(Voulume(2:10,1),PE(2:10,1) ,'-k*');
plot(Voulume(2:10,2),PE(2:10,2) ,'-ko');
plot(Voulume(2:10,3),PE(2:10,3) ,'-ks');
legend('1st order' , '2nd order', '3rd order','Location','NorthWest')
xlabel(' Voulume mm^3');ylabel('pracent error');

%% 

% get the phantom data to fit

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 3;      % Second order is good for up to 5 samples
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

boxSize = repmat(phantomP.rSize,1,nDims);


%%
%% select coil gto fits to the phantom data

% We use this algorithm to
nUseCoils = 4;                         % How many coils to use
c = nchoosek(1:16,nUseCoils);          % All the potential combinations of nCoils
Cor = ones(nUseCoils,size(c,1))*100;   % Initiate the Cor to max
% cLocations = zeros(nUseCoils,size(c,1),nUseCoils);
for kk=1:size(c,1)             % loop over coils combinations
    
    % Correlations between the the measured M0 for each of the coils in
    % this combination.
    A = (corrcoef(phantomP.M0_v(:,c(kk,:)))); 
    
    % Sum the abs correlation after correcting for number of coils and the
    % identity correlations
    Cor(nUseCoils,kk) = sum(sum(abs(triu(A) - eye(nUseCoils,nUseCoils))));
    
    % Keep track of the coil set we test
    % cLocation(nUseCoils,kk,1:nUseCoils)=c(kk,:);
end

% Find the minimum correlation (min abs corr give us the set with corr that
% are closer to zero).  Choose those coils.
[v, minCor] = min(Cor(nUseCoils,:)); 
coils       = c(minCor,:);

%% fit with no regularization


NL   = pdBiLinearFit_lsqSeach(phantomP.M0_v(:,coils),phantomP.pBasis);


%% %% fit with T1  regularization
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
Rmatrix(:,2) = double(1/(phantomP.t1(:)*1000));

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ] = ...
    pdX_valdationLoop_2(lambda,kFold,phantomP.M0_v(:,coils),phantomP.pBasis,Rmatrix,[],[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:))) 

% semilogy(lambda,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('CV error'); 
% X_valdationErr(2,:)./min(X_valdationErr(2,:))

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),phantomP.M0_v(:,coils),phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg));

%%  reshape and scale the PD fits
% this our expextation from the phantom 
PD=ones(size(PD_T1reg));
% PD with T1 reg
PD_T1reg  = reshape(NL_T1reg.PD, boxSize);
%scale     = 1/PD_T1reg(1,1,1);
scale     = mean(PD(:)./PD_T1reg(:));

PD_T1reg  = PD_T1reg.*scale;

% PD with out T1 reg
PD_Noreg  = reshape(NL.PD, boxSize);
scale     = mean(PD(:)./PD_Noreg(:));
%scale     = 1/PD_Noreg(1,1,1);
PD_Noreg  = PD_Noreg.*scale;

% for vizalization let clip the fitted PD to be up to 100% off the real
% value
PD_Noreg(PD_Noreg>2)=2;
PD_Noreg(PD_Noreg<0)=0;

%%  make the figure

mrvNewGraphWin;

hold on
plot((1-PD_Noreg(:))*100,'o' ,'MarkerSize',10,'MarkerFaceColor','b')
plot((1-PD_T1reg(:))*100,'or','MarkerSize',10)

ylabel('Pracent error'); xlabel('Voxels');
axis image; axis square
legend('PD estimate without T1 reg','PD estimate with T1 reg','Location','SouthWest')


mrvNewGraphWin;

plot((1-PD_T1reg(:))*100,'or','MarkerSize',10)

ylabel('Pracent error'); xlabel('Voxels');
axis image; axis square
title('PD estimate with T1 reg')

%  the coefficient of determination (R^2) 
CV1=(calccod(PD_T1reg(:),PD(:))/100).^2
CV2=(calccod(PD_Noreg(:),PD(:))/100).^2
