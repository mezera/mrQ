function  [X_valdationErr,   gEstT, resnorm, FitT]=pdX_valdationLoop( lambda1,kFold,M0,pBasis,R1basis,g0,mask,options)
% Fit M0 for coil gain and PD with PD with different weight (lambda1) for PD
% regularization by T1 (linear relations).
%
%  [X_valdationErr   gEstT, resnorm, FitT]= ...
%     pdX_valdationloop(lambda1, kFold, M0, ...
%                       pBasis,R1basis,g0,mask,options)
%
% Calculates the Cross Validation eror for each regularization wight
%
% AM  & BW VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);

% separate the data to Kfold fit and X-Validation part
[holdX, useX] = getKfooldCVvoxel(nVoxels,kFold);

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
    options = optimset('Display','iter',...
        'MaxFunEvals',Inf,...
        'MaxIter',Inf,...
        'TolFun', 1e-6,...
        'TolX', 1e-10,...
        'Algorithm','levenberg-marquardt');
end

if notDefined('mask');mask=ones(size(M0,1),1);end
T=(unique(mask));T=T(find(T>0));

%% loop over lambda1 regularization wight
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));
% FitT = This is an array of structs
%        (kFold,length(lambda1));
resnorm = zeros(kFold,length(lambda1));
gEstT = zeros(size(pBasis,2),Ncoils,kFold,length(lambda1));

% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        useNow=logical(useX(:,jj));
        
        % Searching on the gain parameters, G.
        [gEstT(:,:,jj,ii), resnorm(ii,jj)] = ...
            lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par,double(M0(useNow,:)),...
            double(pBasis(useNow,:)),  length(find(useNow)), Ncoils, R1basis(useNow,:), lambda1(ii),mask(useNow),T),...
            double(g0),[],[],options);
        
        %  calculate X-Validation error
        holdNow=logical(holdX(:,jj));
        
        %Check if the coil coefficent can explain the hold data
        [FitT(jj,ii)] = ...
            pd_CVtest_voxels(gEstT(:,:,jj,ii) ,double(pBasis(holdNow,:)),double(M0(holdNow,:)));
        
        % Change to sum of squares?
        X_valdationErrKfold(1,jj)= sum(abs( FitT(jj,ii).M0prederr (:)));
        X_valdationErrKfold(2,jj)= sum(  FitT(jj,ii).M0prederr (:).^2);
    end
    
    %sum the  X-Validation error for this lambda1
    X_valdationErr(1,ii)=sum(X_valdationErrKfold(1,:));
    X_valdationErr(2,ii)=sum(X_valdationErrKfold(2,:));
    
    
end

