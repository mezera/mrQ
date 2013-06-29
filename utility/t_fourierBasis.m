%% Create a Fourier Basis
%
% Interacting with Fourier data in MATLAB is often useful.
% This script is part of a tutorial series (see isetbio, where this should
% go ultimately) on manipulating data with Matlab's fft functions.
%
% These basis functions can be used for N-dimensional approximations.
%
% BW Copyright Imageval, LLC 2013

%% Let's work on a 3D case because I am thinking about the for MR data

% I am thinking about X as being in the transform domain
% Let's use a box with 16 samples on a side
nSamples = 16;

% Here are example Fourier coefficients.
% Give them each a try.  One at a time.

% u = 2; v = 2; w = 1;
% u = 1; v = 3; w = 1;
% u = 1; v = 1; w = 3;

X = zeros(nSamples,nSamples,nSamples);

% For visualization, we put the data on a mean background of 1
N = numel(X);
X(1,1,1) = N;     % Make the mean 1
X(u,v,w) = N;     % Here is the coefficient
vBasis = abs(ifftn(X));  % Visualize basis function
showMontage(vBasis);
% mean(vBasis(:))


%% To create the basis functions we could loop like this

X = zeros(nSamples,nSamples,nSamples);
N = numel(X);

% Loop over (u,v,w) values you want in the basis set
%
nFreq = 3;
[U V W] = meshgrid(1:nFreq,1:nFreq,1:nFreq);
nCoef = length(U(:));
bSet = zeros(N,nCoef);
for ii = 1:nCoef
    X = zeros(nSamples,nSamples,nSamples);
    X(U(ii),V(ii),W(ii)) = N;     % Here is the coefficient
    bFunction = abs(ifftn(X));
    % showMontage(bFunction);
    bSet(:,ii) = bFunction(:);
end


% bFunction(:)'*bFunction(:)/length(bFunction(:))

%%
tmp = reshape(bSet(:,10) ,nSamples,nSamples,nSamples); 
showMontage(tmp)


