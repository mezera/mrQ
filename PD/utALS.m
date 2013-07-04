%% Unit testing of alternating linear least squares
%
% AM/BW VISTASOFT Team, 2013

%% If you are in the mrQ directory, run this to set the path
addpath(genpath(fullfile(mrqRootPath)));

%% Run pdPolyPhantomOrder to obtain some reasonable testing parameters

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ.  Other choices don't work, I think (BW).
pOrder   = 1;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
oFlag = true;

% printImages = false;
smoothkernel=[];

% tData contains the various variables for comparing data and polynomial
% approximations. We will turn it into a function before long. Variables
% include M0S_v, pBasis, params, SZ
tData = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, oFlag);
% mrvNewGraphWin; imagesc(tData.pBasis);
% tmp = reshape(tData.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))


% The error to the data can be large, but we don't care.  We just want the
% linear parameters
percentError = 100*tData.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)

%% Initiate params 

coilList = 1:2;

% Linear coil gains estimated from the phantom at one point. 
cGains = tData.params(:,coilList);
pBasis = tData.pBasis;

% Useful parameters for later coding
nCoils = length(coilList);
nVoxels = prod(tData.SZ(1:3));


%% Make simulated M0 data by multiplying PD by coil Gains

% We insist that PD(1) is always 1.
PD = ones(nVoxels,1); 
PD = PD/PD(1);
M0 = diag(PD)*pBasis*cGains;

%% Check that we recover when PD and M0 are known and no noise

% The function begins with M0 and a PD estimate.
% It returns a new PD estimate and a coil gain estimate

PDest = PD;

% For known PD, estimate the coil gains
A = diag(PDest)*pBasis;
cGainsEst = pinv(A)*M0;

% Flesh out the gains into voxel space and then solve for PD
G = pBasis*cGainsEst;
PDest = zeros(size(PD));
for ii=1:nVoxels
    % M0(ii,:)' =  G(ii,:)' * PDest
    PDest(ii) = G(ii,:)' \ M0(ii,:)';
end

% We loop, using PDest as the starting point, and continue until
% PDest stabilizes.

mrvNewGraphWin([],'tall');
subplot(2,1,1)
plot(cGains(:), cGainsEst(:),'o');
identityLine;
subplot(2,1,2)
plot(PD(:), PDest(:),'o');
identityLine;

%% Now do the same, but set up for ridge regression, keeping lambda = 0
% This should be perfect.

lambda1 = 0.0;
PDest = PD;

% We would loop, using PDest as the starting point, and continue until
% PDest stabilizes.

% For known PD, estimate the coil gains
% Minimizing || M0 - A * g || + lambda1 g' * g by ridge regression
A = diag(PDest)*pBasis;
ridgeA = (A'*A + lambda1*eye(size(size(A,2))))\A';
cGainsEst = ridgeA*M0;

PDest = zeros(size(PD));
G = pBasis*cGainsEst;
% M0(ii,:)' =  G(ii,:)' * PDest
for ii=1:nVoxels
    PDest(ii) = G(ii,:)' \ M0(ii,:)';
end

% Force first PD estimate to 1 always
% Compensate coil gains
PDest = PDest/PDest(1);
cGainsEst = cGainsEst*PDest(1);

mrvNewGraphWin([],'tall');
subplot(2,1,1)
plot(cGains(:), cGainsEst(:),'o');
identityLine;
subplot(2,1,2)
plot(PD(:), PDest(:),'o');
identityLine;

%% Now with lambda non-zero, we should deviate just a bit

lambda1 = 0.05;
PDest = PD;

% We would loop, using PDest as the starting point, and continue until
% PDest stabilizes.

% For known PD, estimate the coil gains
% Minimizing || M0 - A * g || + lambda1 g' * g by ridge regression
A = diag(PDest)*pBasis;
ridgeA = (A'*A + lambda1*eye(size(size(A,2))))\A';
cGainsEst = ridgeA*M0;

PDest = zeros(size(PD));
G = pBasis*cGainsEst;
% M0(ii,:)' =  G(ii,:)' * PDest
for ii=1:nVoxels
    PDest(ii) = G(ii,:)' \ M0(ii,:)';
end

% Force first PD estimate to 1 always
% Compensate coil gains
PDest = PDest/PDest(1);
cGainsEst = cGainsEst*PDest(1);

mrvNewGraphWin([],'tall');
subplot(2,1,1)
plot(cGains(:), cGainsEst(:),'o');
identityLine;
subplot(2,1,2)
plot(PD(:), PDest(:),'o');
identityLine;


%% Next, initialize with a the right PD
% But loop a few times and check whether it starts to converge

figH = mrvNewGraphWin([],'tall');
lambda1 = 0.05;
PDest = PD;
lambda2 = 0.1;

A = diag(PDest)*pBasis;
ridgeA = (A'*A + lambda1*eye(size(size(A,2))))\A';

% We would loop, using PDest as the starting point, and continue until
% PDest stabilizes.  At the moment, it makes things worse.
for kk=1:20
    % For known PD, estimate the coil gains
    % Minimizing || M0 - A * g || + lambda1 g' * g by ridge regression

    cGainsEst = ridgeA*M0;
    
    
    % Should we be using ridge regression here, too?
    G = pBasis*cGainsEst;
    % M0(ii,:)' =  G(ii,:)' * PDest
    % g = G(ii,:)'; m = M0(ii,:)';
    % ridgeG = (g'*g + lambda2*eye(size(g,2)))\g'
    % ridgeG*m
    for ii=1:nVoxels
        g = G(ii,:)'; m = M0(ii,:)';
        ridgeG = (g'*g + lambda2*eye(size(g,2)))\g';
        PDest(ii) =  ridgeG*m;
        % Plain regression
        % PDest(ii) = G(ii,:)' \ M0(ii,:)';
    end
    
    % Force first PD estimate to 1 always
    % Compensate coil gains
    PDest = PDest/PDest(1);
    cGainsEst = cGainsEst*PDest(1);
    
    figure(figH);
    subplot(2,1,1)
    plot(cGains(:), cGainsEst(:),'o');
    identityLine;
    subplot(2,1,2)
    plot(PD(:), PDest(:),'o');
    identityLine;
    
end

%%
% The idea is to fix PD and estimate cGains, and then fix cGains and
% estimate PD.  

% How do we solve for PD given cGains?

% 
lambda = 0.5;

% First, assume we know but randomize it
PDest = rand(size(PD));

% Solve for G using ridge regression.  Since M0 = diag(PDest)*G

% Solve for the next PD using ridge regression

% Loop a few times


%%
BL = pdBiLinearFit(tData.M0_v(:,coilList),tData.pBasis,Lambda,maxLoops,sCriterion,[],plotFlag,cGains);

% These are the parameters
% OutPut = pdBiLinearFit(M0_v, pBasis,Lambda,maxLoops, ...
%    sCriterion,PD,plotFlag,TruePar)

PDFinal = reshape(BL.PD,tData.SZ(1:3));

% showMontage(PDFinal)


%% Now let;s simulate different shapes of PD with some noise

% Simulated phantom with all 1's.
% PD = ones(nVoxels,1);
% PD = 'single point';
% PD = 'small region';
PD = 'linear slope';
noiseLevel = 0;
[M0SN, M0S, SNR, PDsim]= simM0(tData.M0S_v(:,coilList),PD,noiseLevel,true);

PDsim = reshape(PDsim,tData.SZ(1:3));
showMontage(PDsim);

BLSim = pdBiLinearFit(M0SN, tData.pBasis,...
    Lambda, maxLoops, sCriterion, [], plotFlag, cGains);
showMontage(reshape(BLSim.PD,tData.SZ(1:3)));

%%
PDFinal = reshape(BLSim.PD,tData.SZ(1:3));
showMontage(PDFinal-PDsim)


%% now lets simulate noise and PD
noiseLevel=5;
PD=[];

[M0SN, M0S,SNR, PDsim]=simM0(tData.M0S_v(:,coilList),PD,noiseLevel,PDtype,plotFlag);
PDsim=reshape(PDsim,tData.SZ(1:3));

BLSim=pdBiLinearFit( M0SN,tData.pBasis,Lambda,maxLoops,sCriterion,[],plotFlag,cGains);
PDFinal = reshape(BLSim.PD,tData.SZ(1:3));
PDsim=reshape(PDsim,tData.SZ(1:3));
showMontage(PDFinal-PDsim)




