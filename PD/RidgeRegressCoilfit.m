function g=RidgeRegressCoilfit(PD,Lamda,M0,pBasis)
%            g=RidgeRegressCoilfit(PD,Lamda,M0,Pbasis)
%
%
%
Ph=diag(PD)*pBasis;

g=(Ph'*Ph+Lamda*eye(size(Ph,2)))^-1 * Ph'*M0;

g=g./g(1);