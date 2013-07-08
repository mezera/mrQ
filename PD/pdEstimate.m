function [PD, G] = pdEstimate(M0_v, pBasis, g)
%pdEstimate - 
%   Detailed explanation goes here
%
% Input 
%  M0_v:    Columns containing the 3D M0 values in a vector
%  pBasis:  Polynomial basis (nPositions x nCoef)
%  g:       Gain coefficients of the polynomials (nCoef x nCoils)
%
% AM/BW Copyright VISTASOFT Team 2013


% This is the coil gains (positions by coils)
G = pBasis*g;

%% Divide the M0_v by the Gains to estimate the PDs
% This is OK, but not perfectly justified because taking the mean is a hack.
% PD = M0_v ./ G;
% PD = mean(PD,2);

%% This is an alternative
% M0_v' = Gn' * diag(PD)
% So, we could solve for PD as a series of linear equations with the first
% column of MO_v'(:,ii) = Gn'(:,ii)*PD(ii)
% So, we loop on ii (locations) 
%   PD(ii) = Gn'(:,ii)\M0_v'(:,ii)
% This could also be a ridge regression, not a pseudoinverse.


%% This is a big matrix approach

% nCoils = size(M0_v,2);
nPositions = size(M0_v,1);

% M0_v has each coil M0 in a column.  The first row is all the coils at
% position 1. The gains of each coils are in the columns of G and the gains
% at position 1 are G(1,:);
% So M0(1,:)' = G(1,:)'*PD(1)
PD = zeros(nPositions,1);
for ii=1:nPositions
    PD(ii) = G(ii,:)' \ M0_v(ii,:)';
end

% Conceptually the same, but much slower
%
% Gdiag = zeros(nPositions*nCoils,nPositions);
% for ii=1:nCoils
%     sRow = (ii-1)*nPositions + 1;
%     eRow = sRow + (nPositions-1);
%     Gdiag(sRow:eRow,:) = diag(G(:,ii));
% end
% % mrvNewGraphWin; imagesc(Gdiag); colormap(gray)
% 
% % This could be a ridge regression also, not a pseudoinverse
% %
% PD = Gdiag \  M0_v(:);


end

