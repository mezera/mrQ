function [PD, G] = pdEstimateT1Reg(M0_v, pBasis, g,Lambda,R1)
% [PD, G] = pdEstimateT1Reg(M0_v, pBasis, g,Lambda,R1)
%
%pdEstimateT1Reg - 
%   Detailed explanation goes here
%
% Input 
%  M0_v:    olumns containing the 3D M0 values in a vector
%  pBasis:  Polynomial basis (nPositions x nCoef)
%  g:           Gain coefficients of the polynomials (nCoef x nCoils)
%Lambda the wight on the T1 linearity regularization
%R1:         1/T1 values (nPositions,1);
%
% AM/BW Copyright VISTASOFT Team 2013

% This is the coil gains (positions by coils)


if notDefined('Lambda');    Lambda=0; end
nPositions = size(M0_v,1);
nCoils = size(M0_v,2);

if notDefined('R1');  
    Lambda=0; 
else
W=eye(nPositions)-R1*R1'/norm(R1);

end

G = pBasis*g;

for ii=1:nCoils
PD(:,ii) = (G(:,ii)'*G(:,ii) + Lambda*W)^-1 * (G(:,ii)'*M0_v(:,ii));
end

%PD = (G'*G + Lambda*W)^-1 * G'*M0_v;

%PD = (G*G' + Lambda*W)^-1 * G*M0_v';


% PD = zeros(nPositions,1);
% for ii=1:nPositions
%     PD(ii) = G(ii,:)' \ M0_v(ii,:)';
%  end