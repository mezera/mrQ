function  [X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ]=pdX_valdationRidgeLoop( lambda1,kFold,M0,pBasis,g0,options,D)
% Fit M0 for coil gain and PD with PD with different weight (lambda1) for
% PD ridge  regularization (minimaized the coefisents)
%
%   [X_valdationErr,   gEstT, resnorm, FitT, useX, kFold ].=..
%  pdX_valdationRidgeLoop( lambda1,kFold,M0,pBasis,g0,options,D)
%
% Calculates the Cross Validation eror for each regularization wight
%
% AM  & BW VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);
nPolyCoef=size(pBasis,2);
% separate the data to Kfold fit and X-Validation part
[ useX, kFold] = getKfooldCVvoxel_full(nVoxels,Ncoils,kFold);

if notDefined('g0')
      PDinit = sqrt(sum(M0.^2,2));    % Sum of squares
      
    % get initial guess
    G  = zeros(nVoxels,Ncoils);
    g0 = zeros(nPolyCoef,Ncoils);
    for ii=1:Ncoils
        G(:,ii)  = M0(:,ii) ./ PDinit(:);         % Raw estimate
        g0(:,ii) = pBasis(:,:) \ G(:,ii);  % Polynomial approximation
    end
    clear G;
end

if notDefined('options')
    options = optimset('Display','off',...  %'iter'final
        'MaxFunEvals',Inf,...
        'MaxIter',200,...
        'TolFun', 1e-6,...
        'TolX', 1e-6,...
        'Algorithm','levenberg-marquardt');
end

%%  Initial the ridge coefigents

if notDefined('D')
      D=ones(nPolyCoef,Ncoils); 
      D(1,:)=0; % don't work on the offset
end

%% loop over lambda1 regularization wight
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));
% FitT = This is an array of structs
%        (kFold,length(lambda1));
resnorm = zeros(kFold,length(lambda1));
gEstT = zeros(nPolyCoef,Ncoils,kFold,length(lambda1));


% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        %select the position to estimate the function
        FitMask=zeros(size(M0));FitMask(find(useX~=jj))=1;FitMask=logical(FitMask);
        
        % Searching on the gain parameters, G.
         W=D*lambda1(ii); %the wights of the regularization
         [gEstT(:,:,jj,ii), resnorm(ii,jj)] = ...
    lsqnonlin(@(par) errFitRidgeNestBiLinear(par,double(M0),double(pBasis),nVoxels,Ncoils,W), ...
             double(g0),[],[],options);
        
        %  calculate X-Validation error
        Xmask=zeros(size(M0));Xmask(find(useX==jj))=1;Xmask=logical(Xmask);
        
        [FitT(jj,ii).err_X, FitT(jj,ii).err_F] = X_validation_errHoldVoxel_full(gEstT(:,:,jj,ii),M0,pBasis,nVoxels,Ncoils,FitMask,Xmask);
        %Check if the coil coefficent can explain the hold data
        
        % Change to sum of squares? not sure which is better
        X_valdationErrKfold(1,jj)= sum(abs( FitT(jj,ii).err_X (:)));
        X_valdationErrKfold(2,jj)= sum(  FitT(jj,ii).err_X (:).^2);
    end
    
    %sum the  X-Validation error for this lambda1
    X_valdationErr(1,ii)=sum(X_valdationErrKfold(1,:));
    X_valdationErr(2,ii)=sum(X_valdationErrKfold(2,:));
    
    
end

