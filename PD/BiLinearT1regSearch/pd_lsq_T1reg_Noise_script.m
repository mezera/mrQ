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
 
Par = OutPut.params(:,[1:5]);

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

noiseLevel = 2;   % ?? Units???
[OutPutSim]= simSPGRs(G,PD,[],[],[],[],noiseLevel,true);

R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = OutPutSim.R1(:);

%% 4)intiate the search parameters

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

%% Loop over  Lambda and check X-validation over voxels
%
% The voxels come from multiple coils.


HoldForCV = 0.2;

% Get X-validation voxels to hold out and to use
[hold use] = getCVvoxel(nVoxels,HoldForCV);

lambda1 = [1e4 1e3 1e2 1e1 1e0 1e-1 0] ;   % Weight on T1 regularization
%% Loop over Lambda 

clist = [1 :5 ];  % Coil list
clear gEst resnorm

for ii=1:length(lambda1),
    % Searching on the gain parameters, G.
    [gEst(:,:,ii), resnorm(ii), dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearT1reg(par, OutPutSim.M0SN(use,clist),...
        OutPut.pBasis(use,:), length(use), length(clist), R1basis(use,:), lambda1(ii)),...
        double(g0(:,clist)),[],[],options);
    
end


%% Test for X-validation error on the held out voxels
close all
plotFlag  = 0;
TruePar   = Par;
TruePD    = OutPutSim.PD(:);
%check for hold voxels
clear Fit
for ii=1:length(lambda1)
    %Check if the coil coefficent can explain the hold data
    Fit{ii} = ...
        pd_CVtest_voxels(gEst(:,:,ii),OutPut.pBasis(hold,:),...
        OutPutSim.M0SN(hold,clist),plotFlag,TruePar(:,clist),TruePD(hold),ii);
end


%% X-validation plot

clear SumR_err Sum_err M_err RMSE
for  ii=1:length(lambda1)

SumR_err(:,ii)=sum(abs(Fit{ii}.M0prederr./OutPutSim.M0SN(hold,clist)));
Sum_err(:,ii)=sum(abs(Fit{ii}.M0prederr));

M_err(:,ii)= Fit{ii}.Meanerr;
RMSE(ii)=Fit{ii}.RMSE;

end

figure;
subplot(1,2,1);
plot(lambda1,sum(Sum_err),'-*')
ylabel('sum of abs err  X-validation');xlabel('lambda')

xlim([-max(lambda1)*0.1 max(lambda1)*1.1])
subplot(1,2,2);
;plot(lambda1,RMSE,'-*')
ylabel('PD RMSE in X-validation voxel');xlabel('lambda')
xlim([-max(lambda1)*0.1 max(lambda1)*1.1])
ylim([-max(RMSE)*0.1 max(RMSE)*1.1])

% figure;plot(lambda1,sum(SumR_err),'-*')
% figure;plot(lambda1,sum(M_err),'-*')


%% visualize the full PD error in space
X_validationScoor=sum(Sum_err)
best=find(X_validationScoor==min(X_validationScoor));
[PDfit, Gfit] = pdEstimate(OutPutSim.M0SN(:,clist), OutPut.pBasis, gEst(:,:,best));
PDfit=PDfit./mean(PDfit(:));
PDsim=OutPutSim.PD;
PDsim=PDsim./mean(PDsim(:));

err=(PDfit-PDsim)./PDsim;

err=reshape(err,OutPut.SZ(1:3));

showMontage(err)

%% the tissue spesipisity to the linearity argumant

%% Loop over Lambda 

clist = [1 :5 ];  % Coil list
clear gEst resnorm

for ii=1:length(lambda1),
    % Searching on the gain parameters, G.
    [gEstT(:,:,ii), resnorm(ii), dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par, OutPutSim.M0SN(use,clist),...
        OutPut.pBasis(use,:), length(use), length(clist), R1basis(use,:), lambda1(ii),OutPutSim.mask(use)),...
        double(g0(:,clist)),[],[],options);
    
end


%% Test for X-validation error on the held out voxels
close all
plotFlag  = 0;
TruePar   = Par;
TruePD    = OutPutSim.PD(:);
%check for hold voxels
clear Fit
for ii=1:length(lambda1)
    %Check if the coil coefficent can explain the hold data
    FitT{ii} = ...
        pd_CVtest_voxels(gEstT(:,:,ii),OutPut.pBasis(hold,:),...
        OutPutSim.M0SN(hold,clist),plotFlag,TruePar(:,clist),TruePD(hold),ii);
end


%% X-validation plot

clear SumR_err Sum_err M_err RMSE
for  ii=1:length(lambda1)

SumR_err(:,ii)=sum(abs(FitT{ii}.M0prederr./OutPutSim.M0SN(hold,clist)));
Sum_err(:,ii)=sum(abs(FitT{ii}.M0prederr));

M_err(:,ii)= FitT{ii}.Meanerr;
RMSE(ii)=FitT{ii}.RMSE;

end

figure;
subplot(1,2,1);
plot(lambda1,sum(Sum_err),'-*')
ylabel('sum of abs err  X-validation');xlabel('lambda')

xlim([-max(lambda1)*0.1 max(lambda1)*1.1])
subplot(1,2,2);
;plot(lambda1,RMSE,'-*')
ylabel('PD RMSE in X-validation voxel');xlabel('lambda')
xlim([-max(lambda1)*0.1 max(lambda1)*1.1])
ylim([-max(RMSE)*0.1 max(RMSE)*1.1])

% figure;plot(lambda1,sum(SumR_err),'-*')
% figure;plot(lambda1,sum(M_err),'-*')


%% visualize the full PD error in space
X_validationScoor=sum(Sum_err)
best=find(X_validationScoor==min(X_validationScoor));
[PDfit, Gfit] = pdEstimate(OutPutSim.M0SN(:,clist), OutPut.pBasis, gEstT(:,:,best));
PDfit=PDfit./mean(PDfit(:));
PDsim=OutPutSim.PD;
PDsim=PDsim./mean(PDsim(:));

err=(PDfit-PDsim)./PDsim;

err=reshape(err,OutPut.SZ(1:3));

showMontage(err)
