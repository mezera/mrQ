function  [X_valdationErr,   gEstT, resnorm, FitT useX, kFold ]=pdX_valdationLoop_2( lambda1,kFold,M0,pBasis,R1basis,g0,mask,options)
% Fit M0 for coil gain and PD with PD with different weight (lambda1) for PD
% regularization by T1 (linear relations).
%
%  [X_valdationErr   gEstT, resnorm, FitT]= ...
%     pdX_valdationloop_2(lambda1, kFold, M0, ...
%                       pBasis,R1basis,g0,mask,options)
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
%[useX ] =getSparceCVvoxel_full(nVoxels,Ncoils);
if notDefined('g0')
    PDinit=Get_PDinit(1,R1basis(:,2));
    %  PDinit = sqrt(sum(M0.^2,2));    % Sum of squares
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

if notDefined('mask');mask=ones(size(M0,1),1);end

%% loop over lambda1 regularization wight
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));
% FitT = This is an array of structs
%        (kFold,length(lambda1));
resnorm = zeros(kFold,length(lambda1));
gEstT = zeros(nPolyCoef,Ncoils,kFold,length(lambda1));

% T is the number of different tissue type. tissue type can be defined by segmtation (WM GM CSF ect...). Each tissue
% type will be used to regularized PD fit by linear model.each tissue type  will have a differnt linear model
T=(unique(mask));T=T(find(T>0));

% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        %select the position to estimate the function
        FitMask=zeros(size(M0));FitMask(find(useX~=jj))=1;FitMask=logical(FitMask);
        
        % Searching on the gain parameters, G.
        [gEstT(:,:,jj,ii), resnorm(ii,jj)] = ...
            lsqnonlin(@(par) errFitNestBiLinearTissueT1reg_full_1(par,double(M0),...
            double(pBasis),  nVoxels, Ncoils, double(R1basis), lambda1(ii),mask,FitMask,T),...
            double(g0),[],[],options);
        
        % the fit err with no Lamda
        %   resnormData(ii,jj)=errFitNestBiLinearTissueT1reg_full(gEstT(:,:,jj,ii),double(M0),...
        %             double(pBasis),  nVoxels, Ncoils, double(R1basis), 0,mask,FitMask);
        
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

