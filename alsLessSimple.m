% Trivial, alternating least squares
% Just like our problem, but smaller.
%
% I have experimented with 1D space. In the case of 1D and 2nd order.  We
% always find a solution that fits the M0. But in this case there are many
% solutions that are close to predicting the data exactly.  We do not
% always get a great PD PDest correlation.  When the correlations are bad,
% we usually have a single voxel that is an outlier.
%
% I have not experimented with the coil gain properties much.
%
% BW

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ.  Other choices don't work, I think (BW).
pOrder   = 2;      % Second order is good for up to 5 samples
nSamples = 5;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
oFlag = true;

% printImages = false;
smoothkernel=[];

% tData contains the various variables for comparing data and polynomial
% approximations. We will turn it into a function before long. Variables
% include M0S_v, pBasis, params, SZ
OutPut = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, [], smoothkernel, oFlag);
% mrvNewGraphWin; imagesc(tData.pBasis);
% tmp = reshape(tData.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))


% The error to the data can be large, but we don't care.  We just want the
% linear parameters
percentError = 100*tData.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)


%% Make small number of samples, two coils and constant gain
cGains = OutPut.params(:,1:3);
pBasis = OutPut.pBasis;
nVoxels = size(pBasis,1);           % Number of voxels

% cGains = [1,2,3 ; 1, -1, 3]';     % Different gains
sCoef = 3;                          % Which coefficient should be scaled to 1
cGains = cGains/cGains(sCoef);          
% Various possible PDs
% PD = (1:nVoxels)/nVoxels; PD = PD(:);
PD = ones(nVoxels,1); PD = PD(:);
% PD = rand(nVoxels,1);

% Here are the M0 data
M0 = diag(PD)*pBasis*cGains;
lambda1 = 0.05;
lambda2 = 0.05;

M0 = M0 + 0.1*randn(size(M0));
PDest = rand(size(PD));      % Random start
% PDest = PD;

% We would loop, using PDest as the starting point, and continue until
% PDest stabilizes.  At the moment, it makes things worse.
% Set up the window
figH = mrvNewGraphWin([],'tall');
for kk = 1:100
    
    % Coil gain parameters are estimated first
    A = diag(PDest)*pBasis;
    ridgeA = (A'*A + lambda1*eye(size(size(A,2))))\A';
    
    % For known PD, estimate the coil gains
    % Minimizing || M0 - A * g || + lambda1 g' * g by ridge regression
    cGainsEst = ridgeA*M0;
    
    % We make sure that the first parameter is always 1
    cGainsEst = cGainsEst/cGainsEst(sCoef);
    
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
        %  PDest(ii) = G(ii,:)' \ M0(ii,:)';
    end
    
    figure(figH);
    subplot(3,1,1), plot(cGains(:), cGainsEst(:),'o');
    identityLine;
    subplot(3,1,2), plot(PD(:), PDest(:),'o');
    identityLine;
    
end

M0pred = diag(PDest)*pBasis*cGainsEst;
subplot(3,1,3), plot(M0(:), M0pred(:),'o');
fprintf('RMSE %f\n',std(M0(:) - M0pred(:)))
identityLine;

c = corrcoef(cGains,cGainsEst);
fprintf('Correlation between CG and CGest: %f\n',c(1,2))

c = corrcoef(PD,PDest);
fprintf('Correlation between PD and PDest: %f\n',c(1,2))

c = corrcoef(M0,M0pred);
fprintf('Correlation between M0 and M0est: %f\n',c(1,2))
