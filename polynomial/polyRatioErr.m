function [CoilCoefErr, PDerr,  estMatrix,  ParMatrix, G, M0, PD,  PDspaceErr]=polyRatioErr(est,params,pBasis)
% Calculate errors from the coil gain estimates
%
%  [CoilCoefErr, PDerr,  estMatrix,  ParMatrix, ...
%      G, M0, PD,  PDspaceErr] = ...
%         polyRatioErr(est,params,pBasis)
%
% Inputs
%  est:     Coil gain estimates
%  params:  Not sure
%  pBasis: polynomial basis functions in columns
%
% Returns
%  CoilCoefErr
%  PDerr
%  estMatrix
%  ParMatrix
%  G
%  M0
%  PD
%  PDspaceErr
%
% AM/BW Copyright Vistasoft team, 2013

% Gain parameter estimates
estMatrix = reshape(est,size(params))';

% Maybe these are the gain parameters from the phantom?
ParMatrix = (params./params(1,1))';

% This would then be some kind of percent error of the gains
CoilCoefErr = (ParMatrix-estMatrix)./ParMatrix;

% This looks like the estimated gains, or equivalently the estimated M0
% when the PD are all equal to one.
G  = pBasis*estMatrix';
M0 = pBasis*ParMatrix';

% Standard deviation of the PD errors, just a number?
PDerr = std(mean(G./M0,2));

if size(ParMatrix,2) == 6, nDimensions = 2; %2ed Order 2D
elseif size(ParMatrix,2) == 10, nDimensions = 3; %3ed order 3D
elseif size(ParMatrix,2) == 4, nDimensions = 3;   % 1st Order 3D
end
nVoxels = size(pBasis,1);

switch nDimensions
    case 2
        dim  = sqrt(nVoxels);
        dims = [dim dim];
    case 3
        dim  = round(nVoxels^(1/3));
        dims = [dim dim dim];
end

% Estimate the PD error, but over coils or space or something? (BW)?
PD = reshape(mean(G./M0,2),dims);

% The standard deviation is the error because we implicitly assume PD=1
% everywhere for the phantom?  Is that right (BW)?
ST = std(G./M0,[],2);

PDspaceErr = reshape(ST,dims);

return
