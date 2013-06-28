%% Script analyzes the polynomial order needed for different size boxes
%
%  We evaluate the percent error in the polynomial approximation for each
%  combination of box size and polynomial order.
%
%  We will use this to justify our choice of polynomial order and box size.
%
%  It appears that nSamples =5 and pOrder = 2 is a good choice that
%  provides a polynomial accuracy to about 1.2 percent precision.
%  We only improve to 0.76 percent when we go to pOrder = 3, which adds 10
%  more parameters.
%
%  N.B.  We only evaluate the error when the M0 is at least 500.  We know
%  that the data are very noisy at levels below this (M0 ranges up to
%  3000). So, this number is good for the coil intensities within 
%
% Copyright Vistasoft Team, 2013

%% If you are in the mrQ directory, run this to set the path
addpath(genpath(fullfile(mrqRootPath)));


%% Set up parameters for N realistic coils
% nCoils   = 32;
% nDims    = 3;
% pOrder   = 2;
% nSamples = 5;      % The box is -nSamples:nSamples
% noiseRange = 500;  % This is the smallest level we consider

%% Get M0 sample data from the coil

%4D
% M0 are the phantom data.  SZ is the size.  meanVal are the mean value
% from each coil
sampleLocation = 3;   % This is
[M0, SZ, meanVal ]= phantomGetData(nSamples,sampleLocation);

% Visualize the box.  The order is each panel is sorted by Z.
% Within each panel there are -nSamples:nSamples points
% showMontage(M0)

% Reshape the M0 data to a 2D image
% This has each coil data in a column
% Each column sweeps out the (x,y,z) values for that coil.
% We think it cycles as x, then, y, then z.  So, (1,1,1), (2,1,1), ... and
% then (1, 2, 1), (2, 2, 1), ...
M0_v = reshape(M0, prod(SZ(1:3)), SZ(4));
mrvNewGraphWin; imagesc(M0_v)

%% This is phantom data and we approximate them by polynomials

% Create the basis functions for the polynomials
[pMatrix,s] = polyCreateMatrix(nSamples,pOrder,nDims);
rSize       = length(s);
nVoxels     = rSize^nDims;
nPolyParams = size(pMatrix,2);

% Get the phantom polynomial coefficients assuming the phantom PD equals
% one.  data = pMatrix * params.  So, pMatrix \ data
params = zeros(nPolyParams , SZ(4));
for ii=1:nCoils
    params(:,ii)= pMatrix\M0_v(:,ii);
end

%% We check whether the approximation is accurate for the box

% M0 prediction as a vector for each coil
M0S_v = zeros(nVoxels,nCoils);
for ii=1:nCoils
    M0S_v(:,ii)= pMatrix*params(:,ii);
end

%  std(  M0S_v(:) - M0_v(:) )
lst = M0_v > noiseRange;
std( (M0S_v(lst) - M0_v(lst)) ./ M0_v(lst))

mrvNewGraphWin; plot(M0S_v(:),M0_v(:),'.');


%% End