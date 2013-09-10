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

%% intilaized X-Validation and fits

% separate the data to Kfold fit and X-Validation part
[ useX, kFold] = getKfooldCVvoxel_full(nVoxels,Ncoils,kFold);

% argumante to save the X-Validation results
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));

% argumante to save the estimated parameters 
resnorm = zeros(kFold,length(lambda1));
gEstT = zeros(nPolyCoef,Ncoils,kFold,length(lambda1));

% T is the number of different tissue type. tissue type can be defined by segmtation (WM GM CSF ect...). Each tissue
% type will be used to regularized PD fit by linear model.each tissue type  will have a differnt linear model
T=(unique(mask));T=T(find(T>0));
%% loop over lambda1 regularization wight
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
        
     
        %%  calculate X-Validation error:
        
        Xmask=zeros(size(M0));Xmask(find(useX==jj))=1;Xmask=logical(Xmask); % Xmask are the hold position form X-Validation.
        
                %Check if the coil coefficent can explain the hold data
        [FitT(jj,ii).err_X, FitT(jj,ii).err_F] = X_validation_errHoldVoxel_full(gEstT(:,:,jj,ii),M0,pBasis,nVoxels,Ncoils,FitMask,Xmask);
        
        % Two posible error function. we don't find a big different between
        % them
        X_valdationErrKfold(1,jj)= sum(abs( FitT(jj,ii).err_X (:))); %sum of absulot erro
        X_valdationErrKfold(2,jj)= sum(  FitT(jj,ii).err_X (:).^2); %RMSE
    end
    
    %sum over the Kfold X-Validation error for this lambda1
    X_valdationErr(1,ii)=sum(X_valdationErrKfold(1,:));
    X_valdationErr(2,ii)=sum(X_valdationErrKfold(2,:));
    
    
end

