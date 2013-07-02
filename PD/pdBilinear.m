%% Illustrate the bilinear alternative fit using ridge regression
%
% AM/BW VISTASOFT Team, 2013

%% If you are in the mrQ directory, run this to set the path
addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 1;      % Second order is good for up to 5 samples
nSamples = 4;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
oFlag = true;

printImages = false;
smoothkernel=[];
% This produces the key variables for comparing data and polynomial
% approximations. We will turn it into a function before long.
% Variables include M0S_v, pBasis, params, SZ
[OutPut] = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, oFlag);
% mrvNewGraphWin; imagesc(OutPut.pBasis);
% tmp = reshape(OutPut.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))

percentError = 100*OutPut.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)

%% Initiate params for the ridge regression
% 
k = 0;% number of iteration
tryagain = 1; % go for the while loop

%  the coil we use. i try different combintion and less combination it is
%  almost as good
coilList = 1:2; 

% the original estimated parameters from the noisy phantom data
Par = OutPut.params(:,coilList);  % Each coil 
Par = Par./Par(1);
        
% Set lambda.
Lambda = 0.1;

% I check lambda's any thing between 1- 0.1 was
% good. smaller lamda convarge slower
% higer lamda also convarge but tolarant some bias. in particular
% when less coils are used.

% The M0 values from the phantom
M0_v = OutPut.M0_v(:,coilList);

%% Let's intialize the PD equal to the sum of sqrs 
% This is the PD that will be returned if we won't have multi coil
% information

PD = sqrt(sum(M0_v.^2,2));
PD = PD ./ mean(PD);   % THere are options - e.g., we could set PD(1) to 1
nCoils    = length(coilList);
nPolyCoef = size(OutPut.pBasis,2);

% When PD is given and M0 is given, we can form an estimate of the coil
% gains.  But note that this estimate is not immediately inside of the
% polynomial space.  So, 
G = zeros(size(M0_v));
g0 = zeros(nPolyCoef,nCoils);
for ii=1:length(coilList)
    G(:,ii)  = M0_v(:,ii) ./ PD;         % Raw estimate
    g0(:,ii) = OutPut.pBasis \ G(:,ii);  % Polynomial approximation
end

% Plot the starting point
mrvNewGraphWin(num2str(k));
subplot(1,2,1); plot(g0(:),Par(:),'.')
identityLine(gca);
xlabel('Estimated gains'); ylabel('True gains');

subplot(1,2,2); plot(PD,'.'); set(gca,'ylim',[.5 1.5]);
ylabel('Estimated PD');xlabel('voxel');

%% Start the loop for solving

k = 0;                % number of iteration
tryagain = 1;         % go for the while loop
maxLoops = 100;
sCriterion = 1e-4;  % Stopping criterion

% M0 is a volume for each coil
% M0_v is a vector for each coil (coils in the columns)

% loop and solve by ridge regration
d = zeros(1,maxLoops);
while tryagain == 1
    k = k+1; % count Ridge regression steps
    
    % fit linear equation coefficients using ridge regression
    g = RidgeRegressCoilfit(PD, Lambda, M0_v, OutPut.pBasis);
    
    % Calculate the coil gain of each coil over space from the estimate
    % Calculate PD for each coil
    [PDn, Gn] = pdEstimate(M0_v, OutPut.pBasis, g);
        
    % normalized to have mean of 1.  Should we also multiply the gains, to
    % keep them consistent?
    PDn = PDn ./ mean(PDn(:));
    Gn  = Gn  .* mean(PDn(:));

    % Check if the new estimate differs from the one before or it's
    % converged
    d(k) = std(PD - PDn);
    if d(k) < sCriterion;
        % If stable to within 1 percent, stop. 
        %
        % We could check the gains, rather than PD, or both
        % if std(G-Gn)<0.01
        
        % if the two solutions are stable --> stop
        tryagain=0;
    else
        % Keep going.
        % Update the new PD and and Gain
        PD = PDn;
        G  = Gn;
        
        % plot the new estimations for coil gain and PD
        mrvNewGraphWin(num2str(k));
        subplot(1,2,1);plot(g(:)/g(1),Par(:),'.')
        identityLine(gca);
        xlabel('Estimated gains'); ylabel('True gains');
        
        subplot(1,2,2);
        plot(PDn,'.'); set(gca,'ylim',[.5 1.5]);
        ylabel('Estimated PD');xlabel('voxel');
        
    end
    % we have to stop some time if it's not convarging
    if k == maxLoops
        tryagain=0;
    end
end

% Average the last two estimates (that look almost the same)
PDFinal = (PD + PDn)/2;
% PDFinal = PDn;

%Reshape PD into a 3D volume and plot
PDFinal = reshape(PDFinal,OutPut.SZ(1:3));
showMontage(PDFinal)
mrvNewGraphWin; semilogy(1:k,d(1:k),'o'); xlabel('steps'); ylabel('Change')
