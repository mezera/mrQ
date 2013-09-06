function OutPut = pdBiLinearFit_lsqSeach(M0,pBasis,PD,options)
% OutPut = pdBiLinearFit_lsqSeach(M0,pBasis,PD)
%
% This fancrtion running an non linear solver for the M0=GPD bi-linear problem.

%
% Input
%  M0:         Columns containing the 3D M0 values in a vector (nVoxels x nCoils)
%  pBasis:       Polynomial basis (nPositions x nCoef)
%  PD            initial PD sulotion
%
%   AM/BW Copyright VISTASOFT Team 2013



if notDefined('options')
options = optimset('Display','iter',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');
end

if notDefined('PD')
    % This is the PD that will be returned if we won't have multi coil
    % information
    PD = sqrt(sum(M0.^2,2));
    PD = PD./ mean(PD);   % There are options - e.g., we could set PD(1) to 1
end
nCoils    = size(M0,2);
nVoxels= size(M0,1);
nPolyCoef = size(pBasis,2);

G  = zeros(nVoxels,nCoils);
g0 = zeros(nPolyCoef,nCoils);

% get initial estimation of the coils
mask1 = ~isnan(PD);   % These are the places we use. 
for ii=1:nCoils
    G(mask1,ii)  = M0(mask1,ii) ./ PD(mask1);         % Raw estimate
    g0(:,ii) = pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end


%% solve the  problen 

[g, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinear(par,M0,pBasis,nVoxels,nCoils)...
         ,double(g0),[],[],options);

    [PD, G] = pdEstimate(M0, pBasis, g);
OutPut.g=g;
OutPut.PD=PD;
OutPut.G=G;
OutPut.resnorm=resnorm;
OutPut.exitflag=exitflag;

