%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% 2) Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 3;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 1;% Which box location
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
%%  initiate the search parameters
clist=[1:5];
nUseCoils=length(clist);
M0=double(OutPut.M0_v(:,clist));
mask=ones(size(M0,1),1);

nVoxels         = size(M0,1);
R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = 1./(OutPut.t1(:)*1000);
R1basis=double(R1basis);
nPolyCoef = size(OutPut.pBasis,2);


options = optimset('Display','iter',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');

%  CHOOSE A START PD
%PDinit = sqrt(sum(M0.^2,2));    % Sum of squares
%  PDinit = rand(size(OutPut.M0_v(:)));   % random
%  PDinit = nan(size(mask)); PDinit(find(mask==1)) = 1; % segmentaion
  PDinit = ones(size(M0,1),1);            %   true solution
PDinit = PDinit(:);
      
% get initial guess
G  = zeros(nVoxels,nUseCoils);
g0 = zeros(nPolyCoef,nUseCoils);

% If the segmentation condition is used, we need to run this. Otherwise,
% there are no NaN values and this doesn't matter.
% We can be specific with what we start the rest will be zeros.
mask1 = ~isnan(PDinit);   % These are the places we use. 
for ii=1:nUseCoils
    G(mask1,ii)  = M0(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) = OutPut.pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end
% lambda1 = [1e4 1e3 1e2 1e1 1e0 1e-1 0] ;   % Weight on T1 regularization
lambda1 = logspace(2,3.5,6);
Kfolod =5;
%% X_valdationLoop
[X_valdationErr ,  X_gEst, Xresnorm, X_Fit]=pdX_valdationLoop( lambda1,Kfolod,M0,  OutPut.pBasis,R1basis,g0,mask,options);


%% Fit all data with the lambda from the best X_valdationErr condition

mrvNewGraphWin;
semilogx(lambda1,X_valdationErr,'*-')
best = find(X_valdationErr==min(X_valdationErr));
fprintf('Choosing lambda(%d) = %f\n',best,lambda(best))

%%
% best = 3;

[gEst, resnorm, dd1, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par,M0,...
    OutPut.pBasis,  nVoxels, length(clist), R1basis, lambda1(best),mask),...
    double(g0),[],[],options);


%%  Visualiztion     
G = OutPut.pBasis*gEst(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
mn = mean(PDfit(:));
PDfit = PDfit/mn;
G = G*mn;

showMontage(PDfit);

% Phantom should have a PD of 1 everywhere
PDsim = ones(size(PDfit));

% Maybe we should aways insist that we have a mean of 1, rather than this.
PDsim = PDsim/mean(PDsim(:)); 

% This is the difference between the simulated and estimated PD, divided by
showMontage( (PDsim - PDfit) ./ PDsim);

% sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))
RMSE = sqrt(mean( (100*(PDsim(:) - PDfit(:))./ PDsim(:) ).^2));
title(['Percent error = ' num2str(RMSE)])

%%
