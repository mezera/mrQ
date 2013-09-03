%% Figure 1
%
% Illustrate how PD and coil gain combine to produce the M0 measurements.
%
% The problem we face is to start with M0 and derive the coil gains.
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


%% Make an image that shows the PD and the PD multiplied by the coil gain
%
% This image shows the effect of coil gain on the M0 data.

%[xx yy]=meshgrid(1:7,1:7)
%PD=xx;

%  Make a PD with some contrast so the image will be visible
%  Could be other things.  This is a sphere with a harmonic and other
%  stuff.  We will finalize later.  It should have different shape compared
%  to the coil gains.  PD is a matrix representing a single slice.
[X,Y] = meshgrid(-nSamples:nSamples,-nSamples:nSamples);
R  = sqrt(X.^2 + Y.^2);
PD = sin(R)./R;PD(isnan(PD))=1;
PD = abs(PD);
PD = sqrt(sqrt(sqrt(PD)));

% Put the coil gains into a matrix. MOS_v is a polynomial fit to a phantom.
% So, we treat these as the coil gain estimates.  The estimates from each
% coil are in the columns of the matrix.  The data in each column
% represents a volume (box) that is -nSamples:nSamples on each side.
G = reshape(phantomP.M0S_v,[phantomP.rSize,phantomP.rSize,phantomP.rSize,nCoils]);

% Create the estimated M0 values in a slice given the PD and coil gains.
M0 = zeros(phantomP.rSize,phantomP.rSize,nCoils);
slice = 3;
for  ii=1:nCoils
    M0(:,:,ii) = PD .* G(:,:,slice,ii);
end

% mrvNewGraphWin;
% for ii=1:nCoils
%     subplot(6,6,ii)
%     imagesc(M0(:,:,ii));
% end
%

% mrvNewGraphWin;
% for ii=1:nCoils
%     subplot(6,6,ii)
%     imagesc(G(:,:,3,ii));
% end

%% Show the PD and a few coil gains.  These combine to make the M0

% PD image
mrvNewGraphWin;
imagesc(PD);
colormap(gray); axis image; axis off
title('PD');
% mrUtilResizeFigure(gcf, 900, 900);
% mrUtilPrintFigure('PD_example_slim.eps');

% Gains
coils = [1 3 5];
n = 1;
mrvNewGraphWin([],'tall');
tmp = G(:,:,slice,coils); mn = min(tmp(:)); mx = max(tmp(:));
for ii=coils
    subplot(3,1,n)
    imagesc(G(:,:,slice,ii));
    caxis([mn mx]);
    colormap(gray); axis image; axis off;
    title(sprintf('Gain for coil %d\n',ii));
    n = n+1;
    
    % mrUtilResizeFigure(gcf, 900, 900);
    % mrUtilPrintFigure(['Gain_example_slice' num2str(ii) '.eps']);
end

% Choose the coils M0
mrvNewGraphWin([],'tall');
n = 1;
tmp = M0(:,:,coils); mn = min(tmp(:)); mx = max(tmp(:));
for ii=coils
    subplot(3,1,n)
    imagesc(M0(:,:,ii)); 
    caxis([mn mx]);
    colormap(gray); axis image; axis off;
    title(sprintf('M0 for coil %d\n',ii));
    n = n+1;
    
    % mrUtilResizeFigure(gcf, 900, 900);
    % mrUtilPrintFigure(['M0_example_slice' num2str(ii) '.eps']);
end

%% End
