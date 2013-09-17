function  [X_valdationErr,   BLFit_RidgeReg, FitT, useX, kFold ]=pdX_valdationLoop_RidgeReg( lambda1,kFold,M0,pBasis,GainPolyPar,maxLoops,sCriterion)
% Fit M0 for coil gain and PD with PD with different weight (lambda1) for PD
% ridge by T1 (linear relations).
%
%  [X_valdationErr   gEstT, resnorm, FitT]= ...
%    pdX_valdationLoop_RidgeReg( lambda1,kFold,M0,pBasis,GainPolyPar,maxLoops,sCriterion)
%
% Calculates the Cross Validation eror for each regularization wight
%
% AM  & BW VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);
nPolyCoef=size(pBasis,2);



%% intilaized X-Validation and fits

% separate the data to Kfold fit and X-Validation part
[holdX useX] =getKfooldCVvoxel(nVoxels,kFold);
% argumante to save the X-Validation results
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));


%% loop over lambda1 regularization wight
% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        %select the position to estimate the function
       
        % Searching on the gain parameters, G.
        mask=logical(holdX(:,jj));
         BLFit_RidgeReg(jj,ii) = pdBiLinearFit(double(M0(mask,:)), pBasis(mask,:), ...
                 1, maxLoops, sCriterion, [], 0 ,GainPolyPar);
      
           
        %%  calculate X-Validation error:
         
                %Check if the coil coefficent can explain the hold data
        [FitT(jj,ii).err_X] = X_validation_errHoldVoxel_Ridge(BLFit_RidgeReg(jj,ii).g ,M0,pBasis,nVoxels,Ncoils,useX(:,jj));
        
        % Two posible error function. we don't find a big different between
        % them
        X_valdationErrKfold(1,jj)= sum(abs( FitT(jj,ii).err_X (:))); %sum of absulot erro
        X_valdationErrKfold(2,jj)= sum(  FitT(jj,ii).err_X (:).^2); %RMSE
    end
    
    %sum over the Kfold X-Validation error for this lambda1
    X_valdationErr(1,ii)=sum(X_valdationErrKfold(1,:));
    X_valdationErr(2,ii)=sum(X_valdationErrKfold(2,:));
    
    
end

