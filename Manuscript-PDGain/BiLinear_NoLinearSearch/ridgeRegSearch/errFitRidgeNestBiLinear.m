function err = errFitRidgeNestBiLinear(g,M0,pBasis,nVoxels,nCoils,W)
%Error function for bilinear coil gain estimate 
%
%  err = errFitRidgeNestBiLinear(g,M0,pBasis,nVoxels,nCoils,,W)
%
% g       - Coil gain polynomial coefficients (nCoeffs x nCoils)
% M0      - Measurement (nVoxels x nCoils)
% pBasis  - Polynomial basis in columns (nVoxels x nCoeffs)
% nVoxels - Number of voxels
% nCoils  - Number of measurement coils
% W       - The ridge regression coefficients
%
% This function is an lsq version of the bilinear ridge regression
% solution. this regularization can be usful to reduce the noie effect.
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

%  the regularization term
reg= W.*sqrt((g).^2);
%The difference between measured and predicted term
err = M0 - M0P ;
% join the two error terms
err=[err(:)'  reg(:)'];

end