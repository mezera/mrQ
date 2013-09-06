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

% This produces the key parameters for the polynomial approximations. We
% will turn it into a function before long. The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);

%% simulte PD
[X,Y Z] = meshgrid(-nSamples:nSamples,-nSamples:nSamples, -nSamples:nSamples);
R  = sqrt(X.^2 + Y.^2 + + Z.^2);
PD = sin(R)./R;PD(isnan(PD))=1;
PD = abs(PD);
PD = sqrt(sqrt(sqrt(PD)));

%% simultae coil Gain (we are using the poylnomyal fits to the phantom data a typical coil function)
%select a set of coils
coils = [1 3 5 9];
% get those coil poylnomyal coeficents 
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multipal the polynomyals coeficents
% and the polynomyal basis.
G = phantomP.pBasis*GainPolyPar;


%% simulte MRI SPGR  signal with and with out Noise
noiseLevel = 2;   % ?? Units???
% simultate the M0 and T1 fits of multi SPGR images. 
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,true);
% MR_Sim is a stracture with multipal fielled that inculde the simulation
% inputs MR sigunalinputsand the calculate of this signal after
% fitting the signal eqation.

%% Try to sparate the M0 component (PD,Gain) using bilinear slosion loop. N0 Noise
 
%set the fiiting loop parmeters 
maxLoops=10e4;
sCriterion = 1e-3;  % Stopping criterion    

% Fit BILinear loop
% the bilinear search is very slow. many many iteration are needed

%     BLFit_N0Noise = pdBiLinearFit(MR_Sim.M0S, phantomP.pBasis, ...
 %  [], maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
%save('/home/avivm/Documents/PD_article/figures/BL_NONoise','BLFit_N0Noise')

%% Alternativly one can use a nonlinear sreach: to reach similar results with very few iteration and much much faster

% Bi-linear problem  non linear solver
NL_noNoise = pdBiLinearFit_lsqSeach(MR_Sim.M0S,phantomP.pBasis);

%% With Noise
% Make PD Gain and M0 with noise
% Fit 
% BLFit_Noise = pdBiLinearFit(MR_Sim.M0SN, phantomP.pBasis, ...
 %   [], maxLoops, sCriterion, [], 1 ,GainPolyPar,PD(:));
% save('/home/avivm/Documents/PD_article/figures/BL_Noise','BLFit_Noise')


%% Alternativly one can use a nonlinear sreach: to reach similar results with very few iteration and much much faster

% Bi-linear problem  non linear solver
NL_Noise = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis);




%% figure  
%% slice of PD the map
slice=4

mrvNewGraphWin
imagesc(PD(:,:,slice));
colormap(gray); axis image; axis off
title('PD');
% mrUtilResizeFigure(gcf, 900, 900);
% mrUtilPrintFigure('PD_example_slim.eps');

%% slice ofPD estimation with out noise form M0
PD_NoNoise=reshape(NL_noNoise.PD,size(PD));
scale=PD(1,1,1)/PD_NoNoise(1,1,1);
PD_NoNoise=PD_NoNoise.*scale;
mrvNewGraphWin
imagesc(PD_NoNoise(:,:,slice));
colormap(gray); axis image; axis off
title('PD estimate with no noise');





%% slice of PD estimation with out noise form M0
PD_Noise=reshape(NL_Noise.PD,size(PD));
scale=PD(1,1,1)/PD_Noise(1,1,1);
PD_Noise=PD_Noise.*scale;
mrvNewGraphWin
imagesc(PD_Noise(:,:,slice));
colormap(gray); axis image; axis off
title('PD estimate with  noise');

%%  PD estimation with out noise vs Simulated PD

MM=minmax([ PD PD_NoNoise]);
mrvNewGraphWin
                hold on
                                                  plot(PD_NoNoise(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')
                xlabel('Estimated PD'); ylabel('True PD');
                                identityLine(gca);
                xlim([MM(1) MM(2)])
                ylim([MM(1) MM(2)])

axis image; axis square  
legend('PD estimate without noise','Location','NorthWest')
                


%%  PD estimation with out noise vs Simulated PD

MM=minmax([PD_Noise PD PD_NoNoise]);
mrvNewGraphWin
                hold on
                                 plot(PD_Noise(:),PD(:),'or','MarkerSize',10)
                                                  plot(PD_NoNoise(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')
                xlabel('Estimated PD'); ylabel('True PD');
                                identityLine(gca);
                xlim([MM(1) MM(2)])
                ylim([MM(1) MM(2)])

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

