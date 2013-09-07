function err = errFitNestBiLinear(g,M0,pBasis,nVoxels,nCoils)
%Error function for bilinear coil gain estimate 
%
%  err = errFitNestBiLinear(g,M0,pBasis,nVoxels,nCoils)
%
% g       - Coil gain polynomial coefficients (nCoeffs x nCoils)
% M0      - Measurement (nVoxels x nCoils)
% pBasis  - Polynomial basis in columns (nVoxels x nCoeffs)
% nVoxels - Number of voxels
% nCoils  - Number of measurement coils
%
% This function works well as an error term for the no-noise case.  When we
% solve in the presence of realistic noise, we use a regularizer.
%
% We know that nVoxels and nCoils are redundant.  We send them in because
% we do not want to call size() every time.  We call this function a lot
% during the search.  Is there a better way to do it?
%
% AM/BW  Vistasoft Team, 2013

%% Estimate coil coeficients
G = pBasis*g;

% Get the best PD for each position a linear sulotion
% this mske  it a nested biliner problem
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0(ii,:)';
end

% Normalize by the first PD value
G  = G  .* PD(1);
PD = PD ./ PD(1);

% Get the predicted M0
M0P = G.*repmat( PD,1,nCoils);

% The difference between measured and predicted
err = M0 - M0P;

end