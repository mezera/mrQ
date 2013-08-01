
% PD run scrip for LSQ search with T1 regolariztion

%% 1) get Poly
%addpath(genpath(fullfile(mrqRootPath)));

%% 2) Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 3;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 5;% Which box location
BasisFlag = 'qr';

printImages = false;
smoothkernel=[];
% This produces the key variables for comparing data and polynomial
% approximations. We will turn it into a function before long.
% Variables include M0S_v, pBasis, params, SZ
[OutPut] = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
% mrvNewGraphWin; imagesc(OutPut.pBasis);
% tmp = reshape(OutPut.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))
percentError = 100*OutPut.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)
%% simulate M0 and R1
 clist=[1:10];
Par = OutPut.params(:,clist);

% Create the coil gains over voxels
G = OutPut.pBasis*Par;
nVoxels         = size(G,1);
nSimulatedCoils = size(G,2);

% Initialize the different spatial structure of the PD
%PD = ones(nVoxels,1);
%PD = 'single point';
%PD = 'small region';
%PD = 'linear slope';
%PD = 'tissue1';
%PD = 'tissue2';  % Subset of voxels
PD = 'phantom';
% Specify the noise level.  Need units

noiseLevel = 1;   % ?? Units???
[OutPutSim]= simSPGRs(G,PD,[],[],[],[],noiseLevel,true);




%% Fit parameters
lambda1 = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0] ;

K=2;
Clist=([3:6]);
pBasis=OutPut.pBasis;
        M0_v =OutPutSim.M0SN(:,:);
        %make a R1 regularazation matrix
        clear R1basis
        R1basis(1:nVoxels,1) = 1; R1basis(:,2) = OutPutSim.R1Fit; R1basis=double(R1basis);
        
        %get the Polinomyal basis
         PDsim = OutPutSim.PD;
         
         
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

        

%% X-validation Fit
  tic
        [X_valdationErrS ,  X_gEstS, XresnormS, X_FitS ]=pdX_valdationLoop_2(lambda1,K,M0_v(:,Clist), pBasis,R1basis,g0); 

 % [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF ]=pdX_valdationLoop_1(lambda1,K,M0_v(:,Clist), pBasis,R1basis); 
%  [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF]=pdX_valdationLoop(lambda1,K,M0_v(:,Clist), pBasis,R1basis,g0); 
            toc    

            best1 = find(X_valdationErrS(1,:)==min(X_valdationErrS(1,:))) % sum of abs err
            best2 = find(X_valdationErrS(2,:)==min(X_valdationErrS(2,:)))% RMSE
%%
   [PDfitS,RMSE1_S]=pdCoilSearch_T1reg( lambda1(best1),M0_v(:,Clist),pBasis,R1basis,X_gEstS(:,:,1,best1),[],[],PDsim);
                if best1~=best2
                    [PDfitS,RMSE2_S]=pdCoilSearch_T1reg( lambda1(best2),M0_v(:,Clist),pBasis,R1basis,X_gEstS(:,:,1,best2),[],[],PDsim);
                else
                    RMSE2_S=RMSE1_S;
                end
%%
 mrvNewGraphWin;
                subplot(1,2,1)

   semilogy(lambda1,X_valdationErrS(1,:),'*-')
                  ylabel('sum of abs err')
                xlabel('lambda')
             title(  ['SIM Choosing lambda is ' num2str(lambda1(best1))  ' gave PDerr ' num2str(RMSE1_S) ])%' for  Ncoils: '  num2str(C)  'and Kfold '  num2str(K) 'for  Nsample:' num2str(nSamples(NS))   ]) 
                subplot(1,2,2)
         semilogy(lambda1,X_valdationErrS(2,:),'*-')
                           title(  ['Choosing lambda is  ' num2str(lambda1(best2))  ' gave PDerr ' num2str(RMSE2_S )]) 

                ylabel('RMSE')
                xlabel('lambda')

 
%% 3D Visualiztion     
PDfitS = reshape(PDfitS,OutPut.SZ(1:3));
showMontage(PDfitS);
title('PD SIM FIT')
PDsim = reshape(PDsim,OutPut.SZ(1:3));

showMontage(PDsim./mean(PDsim(:))-PDfitS./mean(PDfitS(:))  );
title(['Sim fit     RMSE = '   num2str(RMSE1_S) ])



%% now real phantom data
%Clist=([1:4]);

%Clist=([3:6]);

        M0d_v = OutPut.M0_v;
        %make a R1 regularazation matrix

        nVoxels=length(OutPut.t1(:));
        clear R1Dbasis
        R1Dbasis(1:nVoxels,1) = 1; R1Dbasis(:,2) = 1./(OutPut.t1(:)*1000); R1Dbasis=double(R1Dbasis);

        %% visual simlation and real inputs
        mrvNewGraphWin;
        subplot(1,2,1)
        plot(M0_v(:,1),M0d_v(:,1),'*')
        ylabel('M0 coil 1 simlated');xlabel('M0 coil 1 mesuered ')
        identityLine
        subplot(1,2,2)
        plot(OutPut.t1(:)*1000,1./OutPutSim.R1Fit(:),'*')
        ylabel('simlated T1');xlabel('mesuered T1')
        identityLine
        %% Fit parameters
        % all the same beside that we will make a new  initial guess with  the phantom data
        
       % get initial guess
%       nPolyCoef=size(pBasis,2);
%       nCoils=length(Clist);
G  = zeros(nVoxels,nCoils);
g0 = zeros(nPolyCoef,nCoils);
  
        mask1 = ~isnan(PDsim);   % These are the places we use. 
for ii=1:nCoils
    G(mask1,ii)  = M0d_v(mask1,Clist(ii)) ./ PDsim(mask1);         % Raw estimate
    g0(:,ii) = pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end


%% X-validation Fit
  tic
        [X_valdationErrD ,  X_gEstD, XresnormD, X_FitD ]=pdX_valdationLoop_2(lambda1,K,M0d_v(:,Clist), pBasis,R1Dbasis,g0); 

 % [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF ]=pdX_valdationLoop_1(lambda1,K,M0_v(:,Clist), pBasis,R1basis); 
%  [X_valdationErrF ,  X_gEstF, XresnormF, X_FitF]=pdX_valdationLoop(lambda1,K,M0_v(:,Clist), pBasis,R1basis,g0); 
            toc    

            best1 = find(X_valdationErrD(1,:)==min(X_valdationErrD(1,:))) % sum of abs err
            best2 = find(X_valdationErrD(2,:)==min(X_valdationErrD(2,:)))% RMSE
%%
            [PDfitD,RMSE1_D]=pdCoilSearch_T1reg( lambda1(best1),M0d_v(:,Clist),pBasis,R1basis,X_gEstD(:,:,1,best1),[],[],PDsim);
            if best1~=best2
                [PDfitD,RMSE2_D]=pdCoilSearch_T1reg( lambda1(best2),M0d_v(:,Clist),pBasis,R1basis,X_gEstD(:,:,1,best2),[],[],PDsim);
            else
                RMSE2_D=RMSE1_D;
            end
%%
 mrvNewGraphWin;
                subplot(1,2,1)

   semilogy(lambda1,X_valdationErrD(1,:),'*-')
                  ylabel('sum of abs err')
                xlabel('lambda')
             title(  ['DATA Choosing lambda is ' num2str(lambda1(best1))  ' gave PDerr ' num2str(RMSE1_D) ])%' for  Ncoils: '  num2str(C)  'and Kfold '  num2str(K) 'for  Nsample:' num2str(nSamples(NS))   ]) 
                subplot(1,2,2)
         semilogy(lambda1,X_valdationErrD(2,:),'*-')
                           title(  ['Choosing lambda is  ' num2str(lambda1(best2))  ' gave PDerr ' num2str(RMSE2_D )]) 

                ylabel('RMSE')
                xlabel('lambda')

 
%% 3D Visualiztion     
PDfitD = reshape(PDfitD,OutPut.SZ(1:3));
showMontage(PDfitD);
title('PD Phantom FIT')

PDsim = reshape(PDsim,OutPut.SZ(1:3));

showMontage(PDsim./mean(PDsim(:))-PDfitS./mean(PDfitS(:))  );
title(['Phantom    RMSE = '   num2str(RMSE1_D) ])

