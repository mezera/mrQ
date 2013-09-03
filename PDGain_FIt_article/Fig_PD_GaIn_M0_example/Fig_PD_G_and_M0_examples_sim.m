addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
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


%% Sim for figures

%[xx yy]=meshgrid(1:7,1:7)
%PD=xx;

%let's make aPD with shape in space
[X,Y] = meshgrid(-3:3,-3:3)
R = sqrt(x.^2 + Y.^2) 
PD = sin(R)./R;PD(isnan(PD))=1;PD=abs(PD);PD=sqrt(sqrt(sqrt(PD)));



% the Gain in 4D
G=reshape(OutPut.M0S_v,[OutPut.rSize,OutPut.rSize,OutPut.rSize,nCoils]);


for  ii=1:nCoils
    M0(:,:,ii)=PD.*G(:,:,3,ii);
end
figure;
for ii=1:nCoils
    subplot(6,6,ii)
    imagesc(M0(:,:,ii));
end
    
figure;
for ii=1:nCoils
    subplot(6,6,ii)
    imagesc(G(:,:,3,ii));
end

%% make figures

% PD
figure;imagesc(PD);colormap gray
mrUtilResizeFigure(gcf, 900, 900);
mrUtilPrintFigure('PD_example_slim.eps');

% M0
coils=[1 3 5 9]

for ii=coils
    figure;imagesc(M0(:,:,ii));colormap gray;axis off
    
mrUtilResizeFigure(gcf, 900, 900);

mrUtilPrintFigure(['M0_example_slice' num2str(ii) '.eps']);
end

% Gains
for ii=coils
    figure;imagesc(G(:,:,3,ii));colormap gray;axis off
    
mrUtilResizeFigure(gcf, 900, 900);

mrUtilPrintFigure(['Gain_example_slice' num2str(ii) '.eps']);
end