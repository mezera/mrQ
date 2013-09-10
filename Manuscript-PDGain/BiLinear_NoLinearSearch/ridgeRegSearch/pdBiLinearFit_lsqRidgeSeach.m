function solution = pdBiLinearFit_lsqRidgeSeach(M0,pBasis,g0,D,W,options)
% solution = pdBiLinearFit_lsqSeach(M0,pBasis,PD)
%
% This function running an nonlinear solver for the bilinear M0=GPD
% problem.
%
% Input
%  M0:         Columns containing the 3D M0 values in a vector (nVoxels x nCoils)
%  pBasis:       Polynomial basis (nPositions x nCoef)
%  PD            initial PD sulotion
%  D               is the ridge regression withes matrix
% W               a wight  on D
% AM/BW Copyright VISTASOFT Team 2013

%% Parameter initialization
if notDefined('options')
    options = optimset('Display','iter',...
        'MaxFunEvals',Inf,...
        'MaxIter',Inf,...
        'TolFun', 1e-6,...
        'TolX', 1e-10,...
        'Algorithm','levenberg-marquardt');
end

nCoils    = size(M0,2);
nVoxels= size(M0,1);
nPolyCoef = size(pBasis,2);


%% Initial estimation of the coils


if notDefined('g0')
    % get initial guess
    G  = zeros(nVoxels,nCoils);
g0 = zeros(nPolyCoef,nCoils);
mask1 = ~isnan(PDinit);   % These are the places we use.
    for ii=1:Ncoils
        G(:,ii)  = M0(mask1,ii) ./ PDinit(mask1);         % Raw estimate
        g0(:,ii) = pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
    end
end

%%  Initial the ridge coefigents

if notDefined('D')
      D=ones(nPolyCoef,nCoils); 
      D(1,:)=0; % don't work on the offset
end

if notDefined('W')
      D=D;
else
    D=D*W;
end



%% Solve the problen
[g, resnorm,dd1,exitflag] = ...
    lsqnonlin(@(par) errFitRidgeNestBiLinear(par,M0,pBasis,nVoxels,nCoils,D), ...
             double(g0),[],[],options);

[PD, G] = pdEstimate(M0, pBasis, g);
solution.g = g;
solution.PD= PD;
solution.G = G;
solution.resnorm  = resnorm;
solution.exitflag = exitflag;

end
