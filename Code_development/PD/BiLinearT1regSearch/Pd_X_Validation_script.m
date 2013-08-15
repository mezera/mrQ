
% PD run scrip for LSQ search with T1 regolariztion

%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% 2) Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
BasisFlag = 'qr';

printImages = false;
smoothkernel=[];
% This produces the key variables for comparing data and polynomial
% approximations. We will turn it into a function before long.
% Variables include M0S_v, pBasis, params, SZ
[OutPut] = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
% mrvNewGraphWin; imagesc(OutPut.pBasis);
% tmp = reshape(OutPut.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))
percentError = 100*OutPut.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)
%% simulate M0 and R1
 clist=[1:10];
Par = OutPut.params(:,clist);

% Create the coil gains over voxels
G = OutPut.pBasis*Par;
nVoxels         = size(G,1);
nSimulatedCoils = size(G,2);

% Initialize the different spatial structure of the PD
%PD = ones(nVoxels,1);
%PD = 'single point';
%PD = 'small region';
%PD = 'linear slope';
PD = 'tissue1';
%PD = 'tissue2';  % Subset of voxels

% Specify the noise level.  Need units

noiseLevel = 4;   % ?? Units???
[OutPutSim]= simSPGRs(G,PD,[],[],[],[],noiseLevel,true);


%%  intiate the search parameters
R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = OutPutSim.R1(:);

nPolyCoef = size(OutPut.pBasis,2);

options = optimset('Display','iter',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');

%  CHOOSE A START PD
PDinit = sqrt(sum(OutPutSim.M0SN.^2,2));    % Sum of squares
%  PDinit = rand(size(PDsim(:)));   % random
%  PDinit = nan(size(mask)); PDinit(find(mask==1)) = 1; % segmentaion
%  PDinit = PDsim(:);               %   true solution
PDinit = PDinit(:);
      
% get initial guess
G  = zeros(nVoxels,nSimulatedCoils);
g0 = zeros(nPolyCoef,nSimulatedCoils);

% If the segmentation condition is used, we need to run this. Otherwise,
% there are no NaN values and this doesn't matter.
% We can be specific with what we start the rest will be zeros.
mask1 = ~isnan(PDinit);   % These are the places we use. 
for ii=1:nSimulatedCoils
    G(mask1,ii)  = OutPutSim.M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) = OutPut.pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end
lambda1 = [1e4 1e3 1e2 1e1 1e0 1e-1 0] ;   % Weight on T1 regularization
Kfolod =10;
%% X_valdationLoop
[X_valdationErr ,  X_gEst, Xresnorm, X_Fit]=pdX_valdationLoop( lambda1,Kfolod,OutPutSim.M0SN,  OutPut.pBasis,R1basis,g0,OutPutSim.mask,options);


%% fit all data on the best X_valdationErr condition
 figure;plot(lambda1,X_valdationErr,'*-')
best=find(X_valdationErr==min(X_valdationErr));

 

      [gEst, resnorm, dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par,OutPutSim.M0SN,...
      OutPut.pBasis,  nVoxels, length(clist), R1basis, lambda1(best),OutPutSim.mask),...
        double(g0),[],[],options);

 
%%  Visualiztion     
G = OutPut.pBasis*gEst(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ OutPutSim.M0SN(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
showMontage(PDfit);

PDsim=OutPutSim.PD;
PDsim = reshape(PDsim,OutPut.SZ(1:3));

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE) ' the err is : ' num2str( resnorm)] )
