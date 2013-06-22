function [CoilCoefErr PDerr  estMatrix  ParMatrix G M0 PD  PDspaceErr]=polyRatioErr(est,params,pMatrix)


% Gain parameter estimates
estMatrix = reshape(est,size(params))';
ParMatrix=(params./params(1,1))';

CoilCoefErr=(ParMatrix-estMatrix)./ParMatrix;


G=pMatrix*estMatrix';

M0=pMatrix*ParMatrix';

PDerr=std(mean(G./M0,2));

if size(ParMatrix,2)==6
    dim=sqrt(size(pMatrix,1));
dims=[dim dim];
elseif size(ParMatrix,2)==10
    dim=round( (size(pMatrix,1))^(1/3));
dims=[dim dim dim];
end
PD=reshape(mean(G./M0,2),dims);
ST=std(G./M0,[],2);
PDspaceErr=reshape(ST,dims);
