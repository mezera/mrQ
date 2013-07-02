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


% This is the coil gains (columns) across positions (rows)
G = pBasis*g;

% This is the giant matrix of diagonals.  We should turn this into a
% function  G = polyGains2Diag(Gn)

%% Divide the M0_v by the Gains to estimate the PDs
% PD = M0_v ./ G;
% PD = mean(PD,2);

%% This is a big matrix approach

nCoils = size(M0_v,2);
nPositions = size(M0_v,1);

Gdiag = zeros(nPositions*nCoils,nPositions);
for ii=1:nCoils
    sRow = (ii-1)*nPositions + 1;
    eRow = sRow + (nPositions-1);
    Gdiag(sRow:eRow,:) = diag(G(:,ii));
end
% mrvNewGraphWin; imagesc(Gdiag); colormap(gray)

% M0_v(:) = Gdiag * PD
% So, PD
% PD = pinv(Gdiag)*M0_v(:);
% 
PD = Gdiag \  M0_v(:);
% plot(PD)

end

