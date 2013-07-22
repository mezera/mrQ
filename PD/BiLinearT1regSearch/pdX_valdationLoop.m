function  [X_valdationErr,   gEstT, resnorm, FitT]=pdX_valdationLoop( lambda1,Kfolod,M0,pBasis,R1basis,g0,mask,options)
%[X_valdationErr   gEstT, resnorm, FitT]=pdX_valdationloop= lambda1,Kfolod,M0,pBasis,R1basis,g0,mask,options)
% this function fit M0 for coil gain and PD with PD with different wight (lambda1) for PD regularization by T1
% (linear relations). and calculate the Cross Validation eror for each
% regularization wight
%
% AM VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);
% speartate the data to Kfold fit and X-Validation part
[holdX useX] =getKfooldCVvoxel(nVoxels,Kfolod);



%% loop over lambda1 regularization wight
for ii=1:length(lambda1),
    % loop over differnt part of the data Kfold Cross Validation
    for jj=1:Kfolod
useNow=logical(useX(:,jj));

        % Searching on the gain parameters, G.
      [gEstT(:,:,jj,ii), resnorm(ii,jj), dd1, exitflag] = ...
        lsqnonlin(@(par) errFitNestBiLinearTissueT1reg(par,M0(useNow,:),...
       pBasis(useNow,:),  length(find(useNow)), Ncoils, R1basis(useNow,:), lambda1(ii),mask(useNow)),...
        double(g0),[],[],options);
    
%  calculate X-Validation error
holdNow=logical(holdX(:,jj));

    %Check if the coil coefficent can explain the hold data
    FitT(jj,ii) = ...
        pd_CVtest_voxels(gEstT(:,:,jj,ii) ,pBasis(holdNow,:),M0(holdNow,:));

X_valdationErrKfold(jj)= sum(sum(abs(FitT(jj,ii).M0prederr)));
    end
%sum the  X-Validation error for this lambda1
    X_valdationErr(ii)=sum(X_valdationErrKfold);
    
    
end

    