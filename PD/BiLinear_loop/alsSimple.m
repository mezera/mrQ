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

%%
nSamples  = 8;
pOrder    = 2;
sDim      = 1;
basisFlag = 'qr';
pBasis = polyCreateMatrix(nSamples,pOrder,sDim,basisFlag);

%% Make small number of samples, two coils and constant gain
nVoxels = size(pBasis,1);           % Number of voxels

% cGains = [10,2 ; 1, -1]';        % First order
cGains = [100,2,3 ; 100, -1, 3]';   % Second order

% Various possible PDs
PD = (1:nVoxels)/nVoxels; % PD = Shuffle(PD);

% PD = ones(nVoxels,1); PD(:) = PD(:)/PD(1);
% PD = rand(nVoxels,1); PD(:) = PD(:)/PD(1);

% Here are the M0 data
M0 = diag(PD)*pBasis*cGains;

noiseLevel = 0.005;
M0 = M0 + noiseLevel*randn(size(M0));
SNR = 20*log10(mean(M0(:))/noiseLevel);

PDest = ones(size(PD)) * mean(M0(:));
% PDest = rand(size(PD));          % Random start
% PDest = PD;                    % Great start
% Set up the window
figH = mrvNewGraphWin([],'tall');


%%  Alternating least squares solution
lambda1 = 0.0;
lambda2 = 0.0;

% We would loop, using PDest as the starting point, and continue until
% PDest stabilizes.  At the moment, it makes things worse.
for kk=1:100
    
    % Coil gain parameters are estimated first
    A = diag(PDest)*pBasis;
    ridgeA = (A'*A + lambda1*eye(size(size(A,2))))\A';
    
    % For known PD, estimate the coil gains
    % Minimizing || M0 - A * g || + lambda1 g' * g by ridge regression
    cGainsEst = ridgeA*M0;
    

    % Should we be using ridge regression here, too?
    G = pBasis*cGainsEst;
%     M0(ii,:)' =  G(ii,:)' * PDest
%     g = G(ii,:)'; m = M0(ii,:)';
%     ridgeG = (g'*g + lambda2*eye(size(g,2)))\g'
%     ridgeG*m
    for ii=1:nVoxels
        g = G(ii,:)'; m = M0(ii,:)';
        ridgeG = (g'*g + lambda2*eye(size(g,2)))\g';
        PDest(ii) =  ridgeG*m;
        % Plain regression
        %         PDest(ii) = G(ii,:)' \ M0(ii,:)';
    end
    
    % We make sure that the first parameter is always 1
    PDest = PDest / PDest(1);
    cGainsEst = cGainsEst * PDest(1); 
    
    figure(figH);
    subplot(3,1,1), plot(cGains(:), cGainsEst(:),'o');

    xlabel('True gain'); ylabel('Est gain');
    identityLine;
    subplot(3,1,2), plot(PD(:), PDest(:),'o');
    xlabel('True PD'); ylabel('Est PD');
    identityLine;
    
    % Stopping criterion
    M0pred = diag(PDest)*pBasis*cGainsEst;
    if std(M0(:) - M0pred(:)) < 10*noiseLevel
        fprintf('Stopping after %d iterations\n',kk);
        break;
    end

end


subplot(3,1,1)
c = corrcoef(cGains,cGainsEst);
title(sprintf('Correlation: %f\n',c(1,2)))

subplot(3,1,2)
c = corrcoef(PD,PDest);
title(sprintf('Correlation: %f\n',c(1,2)))

subplot(3,1,3), plot(M0(:), M0pred(:),'o');
M0pred = diag(PDest)*pBasis*cGainsEst;
c = corrcoef(M0,M0pred);
identityLine;
xlabel('True M0'); ylabel('Est M0');
title(sprintf('SNR: %.3f (db)',SNR));
fprintf('RMSE %f\n',std(M0(:) - M0pred(:)))

