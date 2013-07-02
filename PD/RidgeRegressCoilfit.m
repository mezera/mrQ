function g = RidgeRegressCoilfit(PD, Lambda, M0, pBasis)
% Solve for gain parameters using ridge regression
% 
%  g = RidgeRegressCoilfit(PD,Lamda,M0,Pbasis)
%
% PD:      Vector of estimated proton densities at nPositions
% Lambda:  Weight for the ridge
% M0:      The coil data, nPositions x nCoils
% pBasis:  Polynomial basis
%
% AM/BW VISTASOFT 2013

Phat = diag(PD)*pBasis;

g = (Phat'*Phat + Lambda*eye(size(Phat,2)))^-1 * Phat'*M0;

end