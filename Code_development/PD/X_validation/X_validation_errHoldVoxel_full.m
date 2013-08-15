function [err_X, err_F] = X_validation_errHoldVoxel_full(g,M0,pBasis,nPositions,nCoils,Fmask,Xmask)

%estimate coil coefficients across the volume
G = pBasis*g;

% Estimate the best PD for each position a linear sulotion
% This makes it a nested bilinear problem
PD = zeros(nPositions,1);
for ii=1:nPositions
    use=Fmask(ii,:);
    PD(ii) = G(ii,use)' \ M0(ii,use)';
end


% Normalize by the first PD value
% Could wait until the end to do this.
% G  = G  .* PD(1);
% PD = PD ./ PD(1);

M0P = G.*repmat( PD,1,nCoils);

% get the M0predicted  error

% err_X = ( M0(Xmask) - M0P(Xmask))./M0(Xmask) ;
% err_F =  (M0(Fmask) - M0P(Fmask))./M0(Fmask) ;
err_X = M0(Xmask) - M0P(Xmask);
err_F =  M0(Fmask) - M0P(Fmask);



