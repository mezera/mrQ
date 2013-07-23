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
% AM VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);

% separate the data to Kfold fit and X-Validation part
[holdX useX] = getKfooldCVvoxel(nVoxels,kFold);


%% loop over lambda1 regularization wight
X_valdationErrKfold = zeros(1,kFold);
X_valdationErr  = zeros(1,length(lambda1));
% FitT = This is an array of structs
%        (kFold,length(lambda1));
resnorm = zeros(kFold,length(lambda1));
% gEstT = zeros()

% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        useNow=logical(useX(:,jj));
        
        % Searching on the gain parameters, G.
        [gEstT(:,:,jj,ii), resnorm(ii,jj)] = ...
            lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par,M0(useNow,:),...
            pBasis(useNow,:),  length(find(useNow)), Ncoils, R1basis(useNow,:), lambda1(ii),mask(useNow)),...
            double(g0),[],[],options);
        
        %  calculate X-Validation error
        holdNow=logical(holdX(:,jj));
        
        %Check if the coil coefficent can explain the hold data
        FitT(jj,ii) = ...
            pd_CVtest_voxels(gEstT(:,:,jj,ii) ,pBasis(holdNow,:),M0(holdNow,:));
        
        % Change to sum of squares?
        X_valdationErrKfold(jj)= sum(abs(FitT(jj,ii).M0prederr(:)));
    end
    
    %sum the  X-Validation error for this lambda1
    X_valdationErr(ii)=sum(X_valdationErrKfold);
    
    
end

