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

%% 3) simulate M0 and R1

Par = OutPut.params(:,[1:3]);

% Create the coil gains over voxels
G = OutPut.pBasis*Par;
nVoxels         = size(G,1);
nSimulatedCoils = size(G,2);

% Initialize the different spatial structure of the PD
%PD = ones(nVoxels,1);
%PD = 'single point';
%PD = 'small region';
%PD = 'linear slope';
%PD = 'tissue1';
PD = 'tissue2';  % Subset of voxels

% Specify the noise level.  Need units
noiseLevel = 5;   % ?? Units???

% Run the simulation
[M0SN, M0S, SNR, PDsim, mask]= simM0(G,PD,noiseLevel,true);

% Create the R1 that we use for regularizing the fit.
% This is the typical linear relationship between R1 (1/T1) and PD
% See also simSPGRs.m
R1  = (2.5./PDsim) - 2.26;
R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = R1(:);

% Put this in the format of a block.
PDsim = reshape(PDsim,OutPut.SZ(1:3));

%% 4)intiate the search parameters

nPolyCoef = size(OutPut.pBasis,2);
options = optimset('Display','iter',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');

%  CHOOSE A START PD
PDinit = sqrt(sum(M0SN.^2,2));    % Sum of squares
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
    G(mask1,ii)  = M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) = OutPut.pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end

%% Loop over  Lambda and check X-validation over voxels
%
% The voxels come from multiple coils.

plotFlag  = 1;
TruePar   = Par;
TruePD    = PDsim(:);
HoldForCV = 0.4;

% Get X-validation voxels to hold out and to use
[hold use] = getCVvoxel(nVoxels,HoldForCV);

lambda1 = [1e2 1e1 1e0 1e-1 0] ;   % Weight on T1 regularization
%% Loop over Lambda 

clist = [1];  % Coil list
clear gEst resnorm

for ii=1:length(lambda1),
    % Searching on the gain parameters, G.
    [gEst(:,:,ii), resnorm(ii), dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearT1reg(par, M0SN(use,clist),...
        OutPut.pBasis(use,:), length(use), length(clist), R1basis(use,:), lambda1(ii)),...
        double(g0(:,clist)),[],[],options);
    
end


%% Test for X-validation error on the held out voxels
plotFlag = 1;
TruePar  = Par;
TruePD   = PDsim(:);

%check for hold voxels
clear Fit
for ii=1:length(lambda1)
    %Check if the coil coefficent can explain the hold data
    Fit{ii} = ...
        pd_CVtest_voxels(gEst(:,:,ii),OutPut.pBasis(hold,:),...
        M0SN(hold,clist),plotFlag,TruePar(:,clist),TruePD(hold),ii);
end

%Conclusion:
%
% The first test was starting from the true solution interesting this
% prefers the middle lambda  1e2 over: [1e4  1e3 1e1 1e0] . In terms of PD
% fit this is not the best but it's very close. And it is a very good
% direction with three coils.
% 
% The best is 1e3 (it's probably proportional to the number of coils.
%
%Checks:
%
% 1 -- the normalization of PD and G is important first it change the best
% lambda PD fit is not so good ect. strange
%
% 2 -- levenberg-marquardt is much much faster the default algoritim and
% gets as good result maybe even sightly better.
%
% 3 -- using a diferent starting point true segmentation random and sum of
%square we got the same solution!!!

%% &c loop over Lambda over given coils
clist = [1 2];
for jj=1:length(lambda1),
% Searching on the gain parameters, G.
[gEst(:,:,jj), resnorm(jj), dd1, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearT1reg(par, M0SN(:,clist),...
    OutPut.pBasis(:,:), nVoxels, length(clist), R1basis, lambda1(jj)),...
    double(g0(:,clist)),[],[],options);

end


%% 7b test for X-validation on other coils

clistT = [1 2 3];

% For hold voxels
for jj=1:length(lambda1),
    %get the predicted PD
    G = OutPut.pBasis*gEst(:,:,jj);
    PD = zeros(nVoxels,1);
    for ii=1:nVoxels
        PD(ii) = G(ii,:)' \ M0SN(ii,clist)';
    end
    PDfit = reshape(PD,OutPut.SZ(1:3));
    %check if the predicted PD explains the other coil
    Fit{jj}=pd_CVtest_coils(PDfit,OutPut.pBasis,M0SN(:,clistT),plotFlag,TruePar(:,clistT),jj);
end

%conclusion: this does not work so good because the X-validation is fit the
%same with all parameters. Therefore it will be hard to use that a way to
%select the right lamda and the right solution.

%% 8a)  let repeat the voxel CV but now fit on a error function that also have the leave 1 out approch with in

clist=[1 2 3];
Fcoils{1}=[1 2 ]  ; Fcoils{2}=[1 3 ];Fcoils{3}=[2 3 ];
Tcoils{1}=[3 ]  ; Tcoils{2}=[2 ];Tcoils{3}=[1 ];
for ii=1:length(lambda1),
    % Searching on the gain parameters, G.
    [gEstL1O(:,:,ii), resnorm(ii), dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearL1OT1Reg(par, M0SN(use,clist),...
        OutPut.pBasis(use,:), length(use),Fcoils,Tcoils , R1basis(use,:), lambda1(ii)),...
        double(g0(:,clist)),[],[],options);
end


%% 8b test for CV on the other voxels
plotFlag=1;
TruePar=Par;
TruePD=PDsim(:);


%for hold voxels
for ii=1:length(lambda1),
%cheack if the coil coefficent can explain the hold data
FitL1O{ii}=pd_CVtest_voxels(gEstL1O(:,:,ii),OutPut.pBasis(hold,:),M0SN(hold,clist),plotFlag,TruePar(:,clist),TruePD(hold),ii);
end
%conclostion:
%very similar to the original sulotion with out the leave 1 out the secound sulotion is the best.




%% 9) more realistic simulation of the SPGR signal and fit (T1,M0)
Par = OutPut.params(:,[1:3]);

G = OutPut.pBasis*Par;
nVoxels = size(G,1);
nSimulatedCoils = size(G,2);

%PD = ones(nVoxels,1);
%PD = 'single point';
%PD = 'small region';
%PD = 'linear slope';
%PD = 'tissue1';
PD = 'tissue2';  % Subset of voxels

noiseLevel = 2;   % ?? Units???
[OutPutSim]= simSPGRs(G,PD,[],[],[],[],noiseLevel,true);

%% 9a intiate fit
options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun', 1e-6,'TolX', 1e-10, 'Algorithm','levenberg-marquardt');
      nPolyCoef = size(OutPut.pBasis,2); 
      
%  START PD
% mean of squr
 PDsosq = sqrt(sum(M0SN.^2,2));
 PDinit=PDsosq;
 PDinit=PDinit(:);

%   random
% PDinit = rand(size(PDsim(:)));
% PDinit=PDinit(:);

%   segmentaion
% PDinit=nan(size(mask));
% PDinit(find(mask==1))=1;
% PDinit=PDinit(:);

%   true solution
%PDinit = OutPutSim.PD(:);     
      
% get inital guess
G  = zeros(nVoxels,nSimulatedCoils);
g0 = zeros(nPolyCoef,nSimulatedCoils);

% we can be specific with what we start the rest will be zeros.
mask1 = ~isnan(PDinit);
for ii=1:nSimulatedCoils
    G(mask1,ii)  = M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) = OutPut.pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end

R1basis(1:nVoxels,1) = 1;
R1basis(:,2) = OutPutSim.R1Fit(:);



%% 7a loop over Lambda for sub set of voxel
clist = [1 2 3];
clear gEst resnorm
for ii=1:length(lambda1),
% Searching on the gain parameters, G.
[gEst(:,:,ii), resnorm(ii), dd1, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearT1reg(par, OutPutSim.M0SN(use,clist),...
    OutPut.pBasis(use,:), length(use), length(clist), R1basis(use,:), lambda1(ii)),...
    double(g0(:,clist)),[],[],options);

end


%% 7b test for CV on the other voxels
clear Fit
plotFlag=1;
TruePar=Par;
TruePD=OutPutSim.PD(:);
%check for hold voxels
for ii=1:length(lambda1),
%cheack if the coil coefficent can explain the hold data
Fit{ii}=pd_CVtest_voxels(gEst(:,:,ii),OutPut.pBasis(hold,:), OutPutSim.M0SN(hold,clist),plotFlag,TruePar(:,clist),TruePD(hold),ii);
end


%%

[PDfit, Gfit] = pdEstimate(OutPutSim.M0SN, OutPut.pBasis, gEst(:,:,3));
PDfit=PDfit./mean(PDfit(:));
PDsim=OutPutSim.PD;
PDsim=PDsim./mean(PDsim(:));

err=(PDfit-PDsim)./PDsim;

err=reshape(err,OutPut.SZ(1:3));

showMontage(err)

%% Extra



%% 5) LSQ fit

% Coil list
clist = [1 2];
lambda1 = 1e4;   % Weight on T1 regularization

% Searching on the gain parameters, G.
[gEst, resnorm, dd1, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearT1reg(par, M0SN(:,clist),...
    OutPut.pBasis, nVoxels, length(clist), R1basis, lambda1),...
    double(g0(:,clist)),[],[],options);

%% 6) Visualiztion     
G = OutPut.pBasis*gEst(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0SN(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE) ' the err is : ' num2str( resnorm)] )
