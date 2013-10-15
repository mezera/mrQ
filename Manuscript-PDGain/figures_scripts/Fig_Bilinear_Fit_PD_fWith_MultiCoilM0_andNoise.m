%% Figure 2
%
% Illustrate howPD can be mesure from M0 while estimate the coil gain using
%  bilinear sulotions
%
% The problem we pace is the effect of noise on this solotions
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

%% simulate PD
[PD, R1] = mrQ_simulate_PD('6',phantomP.nVoxels);


%% Simulate coil gain 
% We use the poylnomial fits to the phantom data a typical coil function

% Select a set of coils
coils = [1 3 5 9];
% Get those coil poylnomyal coeficents
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomials
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;


%% Simultae MRI SPGR signal with and with out noise
noiseLevel = 2;   % ?? Units???

% Simultate the M0 and T1 fits of multi SPGR images.
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,true);

% MR_Sim is a structure with multiple field that includes the simulation
% inputs MR sigunal inputs and the calculate of this signal after fitting
% the signal equation.

%% Separate the M0 component into (PD,Gain) using bilinear ALS
%
% The alternating least squares approach is slow.  We did it for a while,
% and it returns the same answer as the much faster search, below.
%
% Set the fiiting loop parmeters
% maxLoops=10e4;
% sCriterion = 1e-3;  % Stopping criterion
%
% Fit Bilinear loop
% The bilinear search is very slow. many many iteration are needed
%
% BLFit_N0Noise = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
            %     [], maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
% save('/home/avivm/Documents/PD_article/figures/BL_NONoise','BLFit_N0Noise')
%
%% With Noise
% Make PD Gain and M0 with noise
% Fit
% BLFit_Noise = pdBiLinearFit(MR_Sim.M0SN, phantomP.pBasis, ...
%   [], maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
% save('/home/avivm/Documents/PD_article/figures/BL_Noise','BLFit_Noise')

%% The nonlinear search returns the same result with very few iterations

% No simulated noise
NL_noNoise = pdBiLinearFit_lsqSeach(MR_Sim.M0S,phantomP.pBasis);

% With simulated noise
NL_Noise   = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis);

%% Show an example slice of PD the map

mrvNewGraphWin;

slice = 4;
imagesc(PD(:,:,slice));
colormap(gray); axis image; axis off
title('Simulated PD');

% mrUtilResizeFigure(gcf, 900, 900);
% mrUtilPrintFigure('PD_example_slim.eps');

%% Slice of PD estimate in zero noise case 

PD_NoNoise  = reshape(NL_noNoise.PD, boxSize);
scale     = mean(PD(:)./PD_NoNoise(:));
%scale       = PD(1,1,1)/PD_NoNoise(1,1,1);
PD_NoNoise  = PD_NoNoise.*scale;
%%
mrvNewGraphWin
imagesc(PD_NoNoise(:,:,slice));
colormap(gray); axis image; axis off
title('PD estimate with no noise');

%% Slice of PD estimation with noise form M0

PD_Noise = reshape(NL_Noise.PD,boxSize);
scale     = mean(PD(:)./PD_Noise(:));
%scale    = PD(1,1,1)/PD_Noise(1,1,1);
PD_Noise = PD_Noise.*scale;
%%
mrvNewGraphWin
imagesc(PD_Noise(:,:,slice));
colormap(gray); axis image; axis off
title('PD estimate with noise');

%%  Scatter plot of PD estimate with out noise vs Simulated PD

MM = minmax([ PD PD_NoNoise]);
mrvNewGraphWin;
hold on
plot(PD_NoNoise(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca);
xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)])
axis image; axis square
legend('PD estimate without noise','Location','NorthWest')

%%  PD estimation with noise vs Simulated PD

mrvNewGraphWin;

MM = minmax([PD_Noise PD PD_NoNoise]);
hold on
plot(PD_Noise(:),PD(:),'or','MarkerSize',10)
plot(PD_NoNoise(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca);
xlim([MM(1) MM(2)]);ylim([MM(1) MM(2)]);
axis image; axis square
legend('PD estimate with noise','PD estimate without noise','Location','NorthWest')

%% The M0 images
mrvNewGraphWin([],'tall');

mn = min(MR_Sim.M0SN(:)); mx = max(MR_Sim.M0SN(:));
for ii=1:4
    subplot(4,1,ii)
    M0=MR_Sim.M0SN(:,ii);
    M0=reshape(M0,size(PD));
    imagesc(M0(:,:,slice));
    caxis([mn mx]);
    colormap(gray); axis image; axis off;
    title(sprintf('M0 for coil %d\n',ii));
    
    % mrUtilResizeFigure(gcf, 900, 900);
    % mrUtilPrintFigure(['M0_example_slice' num2str(ii) '.eps']);
end

%% End
