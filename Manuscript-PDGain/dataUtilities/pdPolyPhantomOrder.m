function phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, noiseRange,sampleLocation,printImages,smoothkernel, BasisFlag)
% Creates a structure with various parameters needed for estimating gain
%
%  phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, noiseRange,
%                sampleLocation, printImages,smoothkernel)
%
%Inputs:
%  nSamples              number of voxel around the center voxal.
%  sampleLocation        location of central voxel. posibale location are peaked and number 1:5is used to chose betwwn them
%  nCoils                the number of coils
%  nDims                 1, 2 or 3 for the dimention of the problem
%  pOrder                the polyoms order to fit to the problem 1 2 3
%  noiseRange            a value of the M0 under we thing the SNR is too low
%  plotImage             make an image defult false
%  BasisFlag             Force normalization or orthonormal of pBasis
%
%Output:
%  phantomP    a structure that includes these fields
%              M0, M0_v, params, M0S_v, VarEx, SZ, meanVal, pBasis,
%              s, rSize, nVoxels, nVoxels
%
% Copyright Vistasoft Team, 2013

%% Set up parameters for N realistic coils

% Which phantom location
if notDefined('sampleLocation'), sampleLocation = 3;  end
% positions are -nSamples:nSamples
if notDefined('nSamples'),     nSamples = 1;end
if notDefined('nCoils'),       nCoils = 1;  end
if notDefined('nDims'),        nDims = 3;   end  % Spatial dimensions
if notDefined('pOrder'),       pOrder =2;   end
if notDefined('noiseRange'),   noiseRange =5;       end
if notDefined('printImages'),  printImages = false; end
if notDefined('smoothkernel'), smoothkernel = [];
else                           phantomP.smoothkernel=smoothkernel;
end
% if notDefined('oFlag'), oFlag = false; end

%% Get M0 sample data from the coil

% 4D
% M0 are the phantom data.  SZ is the size.  meanVal are the mean value
% from each coil
[M0, SZ, meanVal, t1] = phantomGetData(nSamples,sampleLocation,smoothkernel);

% Visualize the box.  The order is each panel is sorted by Z.
% Within each panel there are -nSamples:nSamples points
% showMontage(M0)

% Reshape the M0 data to a 2D image
% This has each coil data in a column
% Each column sweeps out the (x,y,z) values for that coil.
% We think it cycles as x, then, y, then z.  So, (1,1,1), (2,1,1), ... and
% then (1, 2, 1), (2, 2, 1), ...
M0_v = reshape(M0, prod(SZ(1:3)), nCoils);
if printImages == 1
    mrvNewGraphWin; imagesc(M0_v)
end
%% This is phantom data and we approximate them by polynomials

% Create the basis functions for the polynomials.  Orthogoanlized.
[pBasis, s, pTerms, W]  = polyCreateMatrix(nSamples,pOrder,nDims,BasisFlag);
rSize       = length(s);
nVoxels     = rSize^nDims;
nPolyParams = size(pBasis,2);

% Get the phantom polynomial coefficients assuming the phantom PD equals
% one.  data = pBasis * params.  So, pBasis \ data
params = zeros(nPolyParams , nCoils);
for ii=1:nCoils
    params(:,ii)= pBasis \ M0_v(:,ii);
end

%% We check whether the approximation is accurate for the box

% M0 prediction as a vector for each coil
M0S_v = zeros(nVoxels,nCoils);
for ii=1:nCoils
    M0S_v(:,ii)= pBasis*params(:,ii);
end

%  std(  M0S_v(:) - M0_v(:) )
lst = M0_v > noiseRange;
percentError=std( (M0S_v(lst) - M0_v(lst)) ./ M0_v(lst));
if printImages==1
    mrvNewGraphWin; plot(M0S_v(:),M0_v(:),'.');
end

%% Assign parameters to output
phantomP.M0     = M0;
phantomP.M0_v   = M0_v;
phantomP.params = params;
phantomP.M0S_v = M0S_v;
phantomP.percentError = percentError;
phantomP.meanVal= meanVal;
phantomP.pBasis = pBasis;
phantomP.s = s;
phantomP.rSize   = rSize;
phantomP.nVoxels = nVoxels;
phantomP.pTerms = pTerms;
phantomP.SZ  = SZ;
phantomP.W   = W;
phantomP.t1  = t1;

end

