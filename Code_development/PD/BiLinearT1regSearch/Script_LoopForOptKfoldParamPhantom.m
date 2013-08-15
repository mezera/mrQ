%Fit parametes
nDims=3;
BasisFlag = 'qr';
pOrder=3;
% Serach parameters
sampleLocation=[2 3 4 5 1 ];
nSamples=[3 4 5 ];
LastCoil=[2 3 4 5 ];
lambda1 = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0] ;

X_valdationErrAll=zeros(2,length(lambda1),length(sampleLocation),length(nSamples),LastCoil(end),LastCoil(end));
RMSEAll=zeros(2,length(sampleLocation),length(nSamples),LastCoil(end),LastCoil(end));
RunTime=zeros(length(sampleLocation),length(nSamples),LastCoil(end),LastCoil(end));
%%
for Loc=1:sampleLocation
    
    for NS=1:length(nSamples)
        % get M0 and t1
        [M0, SZ, ~, t1]= phantomGetData(nSamples(NS),sampleLocation(Loc));
        M0_v = reshape(M0, prod(SZ(1:3)), SZ(4));
        %make a R1 regularazation matrix
        nVoxels=length(t1(:));
        clear R1basis
        R1basis(1:nVoxels,1) = 1; R1basis(:,2) = 1./(t1(:)*1000); R1basis=double(R1basis);
        %get the Polinomyal basis
        [pBasis, s, pTerms ]  = polyCreateMatrix(nSamples(NS),pOrder,nDims,BasisFlag);
        PDsim = ones(size(M0_v,1),1);% % Phantom should have a PD of 1 everywhere

        for C=LastCoil
            MaxkFold=C;
            for K=2:MaxkFold
                tic;
                    % X validation fit with T1 regularization 
                [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF]=pdX_valdationLoop_1(lambda1,K,M0_v(:,1:C), pBasis,R1basis); 
                
                % minimum X_valdation error
                best1 = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:))); % sum of abs err
                best2 = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)));% RMSE

               %fit with all the data given the best X_valdation and calculate the RMSE
               RunTime(Loc,NS,C,K)=toc;
               
                [PDfit,RMSE1]=pdCoilSearch_T1reg( lambda1(best1),M0_v(:,1:C),pBasis,R1basis,[],[],[],PDsim);
                if best1~=best2
                    [PDfit,RMSE2]=pdCoilSearch_T1reg( lambda1(best2),M0_v(:,1:C),pBasis,R1basis,[],[],[],PDsim);
                else
                    RMSE2=RMSE1;
                end
        
        %save the result         
                X_valdationErrAll(:,:,Loc,NS,C,K)=X_valdationErrF;
                RMSEAll(1,Loc,NS,C,K)=RMSE1;
                RMSEAll(2,Loc,NS,C,K)=RMSE2;
                ['saveing loc: ' num2str(Loc) ' Nsample ' num2str(NS) ' last coil: ' num2str(C) ' kfold: '  num2str(K) ]
 save('X_validationFit_1','RMSEAll','RunTime','X_valdationErrAll','sampleLocation','nSamples','LastCoil','lambda1');
            end
        end
    end
end