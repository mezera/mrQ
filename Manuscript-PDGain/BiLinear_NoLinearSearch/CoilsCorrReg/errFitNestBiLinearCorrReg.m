function err = errFitNestBiLinearCorrReg(g,M0,pBasis,nCoils,RegWeight,TissueMask,coefdat,M0mask,nPositions)
% Bilinear estimation subject to regularization by coil correlations
%
%  err =
%    errFitNestBiLinearCorrReg(g,M0,pBasis,nCoils,RegWeight,TissueMask,nPositions)
%
% AM/BW(c) VISTASOFT Team, 2013

%% 
%estimate coil coefficients across the volume
G = pBasis*g;

% Estimate the best PD for each position a linear sulotion
% This makes it a nested bilinear problem

% PD = zeros(nPositions,1);
%   for ii=1:nPositions
% % %     %use=M0mask(ii,:);
%       PD1(ii) = G(ii,M0mask(ii,:))' \ M0(ii,M0mask(ii,:))';
%   end

PD=M0./G;
PD(~M0mask)=nan;
PD=nanmean(PD,2);

%%
%turn your attention to the correlation (overlap) between the coil
% gain functions.  We compute the corrcoefs and and take out the lower
% triangular (non-redundant) part of these.  -1 means take everything below
% the diagonal
coefG =  tril(corrcoef(G),-1);
% Outside of this routine we already calculated
%   coefdat =tril(corrcoef(M0),-1);
% Here we compare the two correlation coefficients
% The correlation coefficients of the data (coefdat) should be larger than
% the corr coef of the gain (coefG).

corrErr = (coefdat(coefdat~=0)-coefG(coefdat~=0))./abs(coefdat(coefdat~=0)); 
% If the difference is positive, there is no regularization penalty.
corrErr(corrErr>0)=0;

%%
% Normalize by the first PD value
G  = G  .* mean(PD(TissueMask>0));
PD = PD ./ mean(PD(TissueMask>0));

% get the predicted M0 for all of the coils
M0P = G.*repmat( PD,1,nCoils);

%% error
% The error is PD - PDpred
%

% The error is a vector with positive and negative values representing the
% M0 difference and the coil corrlations  failures
err = [ M0(:) - M0P(:); (RegWeight)*(corrErr)];

end