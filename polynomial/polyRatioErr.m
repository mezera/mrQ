function Res = polyRatioErr(est,params,pBasis)
% Calculate errors from the coil gain estimates
%
%  Res = polyRatioErr(est,params,pBasis)
%
% Inputs
%  est:     Coil gain estimates
%  params:  Not sure
%  pBasis: polynomial basis functions in columns
%
% Returns Res, a structure with slots
%   estGainParams  - Estimated polynomial gain coefficients
%   trueGainParams - True polynomial gain coefficients
%   spaceGains     - Gains at each position
%   M0        - Simulated M0
%   PD        - M0 divided by gain
%   PDspaceErr - Proton density differences across coils
%   PDerr     - Proton density error in the estimate
%
% AM/BW Copyright Vistasoft team, 2013

% Gain parameter estimates
Res.estGainParams = reshape(est,size(params))';

% Maybe these are the gain parameters from the phantom?
Res.trueGainParams = (params./params(1,1))';

% This would then be some kind of percent error of the gains
Res.coilCoefErr = (Res.trueGainParams - Res.estGainParams)./ Res.trueGainParams;

% This looks like the estimated gains, or equivalently the estimated M0
% when the PD are all equal to one.
Res.spaceGains  = pBasis*Res.estGainParams';
Res.M0 = pBasis*Res.trueGainParams';

% Standard deviation of the PD errors, just a number?
Res.PDerr = std(mean(Res.spaceGains ./ Res.M0,2));

if size(Res.trueGainParams,2) == 6, nDimensions = 2; %2ed Order 2D
elseif size(Res.trueGainParams,2) == 10, nDimensions = 3; %3ed order 3D
elseif size(Res.trueGainParams,2) == 4, nDimensions = 3;   % 1st Order 3D
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
Res.PD = reshape(mean(Res.spaceGains ./ Res.M0,2),dims);

% The standard deviation is the error because we implicitly assume PD=1
% everywhere for the phantom?  Is that right (BW)?
Res.ST = std(Res.spaceGains ./ Res.M0,[],2);

Res.PDspaceErr = reshape(Res.ST,dims);

return
