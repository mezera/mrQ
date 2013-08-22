addpath(genpath(fullfile(mrqRootPath)));
%%

%Fit parametes
nDims=3;
BasisFlag = 'qr';
pOrder=3;
% Serach parameters
sampleLocation=[2 3 4 5 1 ];
nSamples=[3 4 5 ];
LastCoil=[2 3 4 5 ];
lambda1 = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0] ;
%lambda1 = [ 1e3 0] ;

%smoothkernel=3;

Loc=5;NS=3;C2=3; C1=2;K=2
%%
Clist=([3:6]);

    [M0, SZ, ~, t1]= phantomGetData(nSamples(NS),sampleLocation(Loc));
        M0_v = reshape(M0, prod(SZ(1:3)), SZ(4));
        %make a R1 regularazation matrix
        nVoxels=length(t1(:));
        clear R1basis
        R1basis(1:nVoxels,1) = 1; R1basis(:,2) = 1./(t1(:)*1000); R1basis=double(R1basis);
        %get the Polinomyal basis
        [pBasis, s, pTerms ]  = polyCreateMatrix(nSamples(NS),pOrder,nDims,BasisFlag);
         PDsim = ones(size(M0_v,1),1);% % Phantom should have a PD of 1 everywhere
      % get initial guess
      nPolyCoef=size(pBasis,2);
      nCoils=length(Clist);
G  = zeros(nVoxels,nCoils);
g0 = zeros(nPolyCoef,nCoils);
  
        mask1 = ~isnan(PDsim);   % These are the places we use. 
for ii=1:nCoils
    G(mask1,ii)  = M0_v(mask1,Clist(ii)) ./ PDsim(mask1);         % Raw estimate
    g0(:,ii) = pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end

        
        %%
        
        tic
        [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF ]=pdX_valdationLoop_2(lambda1,K,M0_v(:,Clist), pBasis,R1basis,g0); 

 % [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF ]=pdX_valdationLoop_1(lambda1,K,M0_v(:,Clist), pBasis,R1basis); 
%  [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF]=pdX_valdationLoop(lambda1,K,M0_v(:,Clist), pBasis,R1basis,g0); 
            toc    
                % minimum X_valdation error
                best1 = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:))) % sum of abs err
                best2 = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)))% RMSE

        %%        
                  mrvNewGraphWin;
                subplot(1,2,1)

                 semilogy(lambda1,X_valdationErrF(1,:),'*-')
                best = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:)));
                  ylabel('sum of abs err')
                xlabel('lambda')
              title(  ['Choosing lambda is ' num2str(lambda1(best))  ' gave PDerr ' num2str(RMSE1) ])%' for  Ncoils: '  num2str(C)  'and Kfold '  num2str(K) 'for  Nsample:' num2str(nSamples(NS))   ]) 
                subplot(1,2,2)
         semilogy(lambda1,X_valdationErrF(2,:),'*-')
                best = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)));
                              title(  ['Choosing lambda is  ' num2str(lambda1(best))  ' gave PDerr ' num2str(RMSE2 )]) 

                ylabel('RMSE')
                xlabel('lambda')
                
                %%
                [PDfit,RMSE1]=pdCoilSearch_T1reg( lambda1(best1),M0_v(:,Clist),pBasis,R1basis,X_gEstF(:,:,1,best1),[],[],PDsim);
                if best1~=best2
                    [PDfit,RMSE2]=pdCoilSearch_T1reg( lambda1(best2),M0_v(:,Clist),pBasis,R1basis,X_gEstF(:,:,1,best2),[],[],PDsim);
                else
                    RMSE2=RMSE1;
                end
                                %%

                PDfit=reshape(PDfit(:),SZ(1:3));
                 showMontage(PDfit)
    