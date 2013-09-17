function [PDfit,RMSE,G,gEst,resnorm,exitflag]= ...
     pdCoilSearch_T1reg(lambda1,M0,pBasis,Rmatrix,g0,mask,options,PDsim)
% Search for the coil gain functions using T1 regularization
%
% Inputs:
%  lambda1 - Regularization weight
%  M0      - data
%  pBasis  - poly basis functions
%  Rmatrix - Matrix for solving linear coefficients of PD-R1 relationship
%  g0      - The initial guess on the polynomial coefficients
%  mask    - The location of the data and an integer defining the data
%            types
%
% Returns 
%  PDfit
%  RMSE
%  gEst
%  resnorm
%  exitflag
%
%
% AM Vistasoft Team 2013


%%
nVoxels=size(M0,1);
Ncoils=size(M0,2);
nPolyCoef=size(pBasis,2);
G  = zeros(nVoxels,Ncoils);

if notDefined('g0')
    PDinit = sqrt(sum(M0.^2,2));    % Sum of squares
    % get initial guess
    g0 = zeros(nPolyCoef,Ncoils);
    for ii=1:Ncoils
        G(:,ii)  = M0(:,ii) ./ PDinit(:);         % Raw estimate
        g0(:,ii) = pBasis(:,:) \ G(:,ii);  % Polynomial approximation
    end
end

if notDefined('options')
    options = optimset('Display','iter',...
        'MaxFunEvals',Inf,...
        'MaxIter',Inf,...
        'TolFun', 1e-6,...
        'TolX', 1e-10,...
        'Algorithm','levenberg-marquardt');
end

if notDefined('mask');mask=ones(size(M0,1),1);end

% The T can have several distinct values, referring to different tissue
% types
T = (unique(mask));
T = T(T>0);

%% Fit with T1 regularization

% The error function segments the data segmenting the data in the volume
% according to the tissue type specified in T.  Each type is allowed to
% have its own linear relationship.
[gEst, resnorm, ~, exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par, double(M0),...
     double(pBasis),  nVoxels, Ncoils, double(Rmatrix), lambda1,mask,T),...
    double(g0),[],[],options);


%% Estimate PD and error
G = pBasis*gEst(:,:);
PDfit = zeros(nVoxels,1);
for ii=1:nVoxels
    PDfit(ii) = G(ii,:)' \ M0(ii,:)';
end

% We  insist that we have a mean of 1
mn = mean(PDfit(:));
PDfit = PDfit/mn;
G = G*mn;

if notDefined('PDsim')
    RMSE=[];
else
    %  we  insist that we have a mean of 1
    PDsim = PDsim/mean(PDsim(:));
    RMSE = sqrt(mean( (100*(PDsim(:) - PDfit(:))./ PDsim(:) ).^2));
end


