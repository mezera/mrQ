%%
%conclusions:
%with noise no matter were the starting point is even the true we ended up
%in wrong sulotion in PD but better sulotion in M0. 
%Fit : PD error is 0.21 M0 error 3.23  
% True PD error is 0.0042 M0 error 3.2571
% cheack this sulotion with no noise
%Fit : PD error is 0.21 M0 error 0.4  
%True PD error is 0 M0 error 0
%
% so we are totlay dominated by the noise!!! in error function

%%
% i like to answer 
% 1) if we fit part of the date how it fits the other part?
% 2) if we fit on one set of coils how it predict a differnt set?
%%




%% 1) get Poly
addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
nSamples = 3;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
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

%% 2) simulte M0
Par=OutPut.params(:,[1 :4]);
%Par(1,:)=Par(1,:)./100; % what if we keep the constant close to the other values 
G=OutPut.pBasis*Par;
nVoxels=size(G,1);
nCoilsS=size(G,2);

%PD = ones(nVoxels,1);
% PD = 'single point';
% PD = 'small region';
%PD = 'linear slope';
%PD = 'tissue1';
PD = 'tissue2';

noiseLevel = 5;
[M0SN, M0S, SNR, PDsim, mask]= simM0(G,PD,noiseLevel,true);

PDsim = reshape(PDsim,OutPut.SZ(1:3));

%% intiate the search 
options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun', 1e-6,'TolX', 1e-10);
      nPolyCoef = size(OutPut.pBasis,2); 
           PDsosq = sqrt(sum(M0SN.^2,2));

 PDinit=PDsim(:);
         
         % get inital guess
 G = zeros(nVoxels,nCoilsS);    
g0 = zeros(nPolyCoef,nCoilsS);
% we can be spesipic with what we start the rest will be zeros.
mask1=~isnan(PDinit);
for ii=1:nCoilsS
    G(mask1,ii)  = M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) =OutPut. pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end
  %%  lets fit on different combination is more is better
  clist=[3 4];

[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinear(par,M0SN(:,clist),OutPut.pBasis,nVoxels,length(clist))...
         ,double(g0(:,clist)),[],[],options);
     
G = OutPut.pBasis*res1(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0SN(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE) ' the err is : ' num2str( resnorm)] )


% coil 1 2 got RMSE 0.1681
% coil 1 3 got RMSE   1.2975
% coil 1 4 got RMSE  0.2130
% coil 2 3 got RMSE 0.6873
% coil 2 4 got RMSE 0.1974
% coil 3 4 got RMSE   0.1364

% coil 2 3 4 got RMSE 0.1872
% coil 1 3 4 got RMSE  0.0577
% coil 1 2 3 got RMSE 0.0501
% coil 1 2 4 got RMSE  0.125

%coil 1 2 3 4  got RMSE 0.0681
%% lets fit on two and cheack on the other two
  clist=[3 4];
clist2=[2 1];
[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinear(par,M0SN(:,clist),OutPut.pBasis,nVoxels,length(clist))...
         ,double(g0(:,clist)),[],[],options);
  
G = OutPut.pBasis*res1(:,:);
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0SN(ii,clist)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));


BLSim = pdBiLinearFit_1(M0SN(:,clist2), OutPut.pBasis, ...
    0, 1 ,0, PDfit(:), 1, Par(:,clist2));
 M0Fit=BLSim.M0Fit

BLSim = pdBiLinearFit_1(M0SN(:,clist2), OutPut.pBasis, ...
    0, 1 ,0, PDsim(:), 1, Par(:,clist2));

M0Fittrue=BLSim.M0Fit

%
% the true PD better explain  M0 of coil 3 and 4 then the wrong PD fitted
% if coil 1,2 
%fit on 12 test on 34 TrueErr=3.5329 wrongErr=  3.6260
%fit on 14 test on 32 TrueErr= 3.7330 wrongErr=  3.8022
%fit on 13 test on 42 TrueErr=  3.7526        wrongErr=   6.3258
%fit on 23 test on 41 TrueErr=  3.4921   wrongErr=   3.7648
%fit on 24 test on 13 TrueErr=  3.5199   wrongErr=4.1016
%fit on 34 test on 12 TrueErr= 3.5097   wrongErr= 3.9573
% concultion we fitting noise but maybe we can use the other data to
% control it? but how?
%% let fit on partial set
voxelLocations=randperm(nVoxels);
cut=round(nVoxels*0.5);

CVvox=voxelLocations(1:cut);
Fitvox=voxelLocations(cut+1:end);
%nvoxelf=length(Fitvox);

     
       clist=[1 3 ];

[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinear(par,M0SN(Fitvox,clist),OutPut.pBasis(Fitvox,:),length(Fitvox),length(clist))...
         ,double(g0(:,clist)),[],[],options);
  
     G = OutPut.pBasis(CVvox,:)*res1(:,:);
     NV=length(CVvox);
PDfit = zeros(NV,1);
for ii=1:NV
    PDfit(ii) = G(ii,:)' \ M0SN(CVvox(ii),clist)';
end


BLSim = pdBiLinearFit_1(M0SN(CVvox,clist), OutPut.pBasis(CVvox,:), ...
    0, 1 ,0, PDfit(:), 1, Par(:,clist));
 M0Fit=BLSim.M0Fit

 PDsimCV=PDsim(CVvox);
BLSim = pdBiLinearFit_1(M0SN(CVvox,clist), OutPut.pBasis(CVvox,:), ...
    0, 1 ,0, PDsimCV(:), 1, Par(:,clist));
M0Fittrue=BLSim.M0Fit


%
%
% the true PD better explain  M0 of voxel that was not fit  then the wrong PD fitted
% to other voxels2 
%

% coil 1 2 TrueErr=  3.3314 wrongErr=   3.3942
% coil 1 3 TrueErr= 3.4874   wrongErr= 3.5883
% coil 1 4 TrueErr=  3.3534   wrongErr=  3.5815
% coil 2 3 TrueErr= 3.3534  wrongErr= 3.5937
% coil 2 4 TrueErr= 3.2000   wrongErr=   3.4168
% coil 3 4 TrueErr=3.3975   wrongErr= 3.6776
%
% coil 2 3 4 TrueErr=4.1907   wrongErr=  4.2541
% coil 1 3 4 TrueErr= 4.0248   wrongErr= 4.0703
% coil 1 2 3 TrueErr= 4.0928   wrongErr=  4.1680
% coil 1 2 4 TrueErr=4.0063   wrongErr= 4.0028   --> this one go the oposit
% way!!! (the are besiclly the same) secound run random cut again close but
% oposite   TrueErr=4.0496   wrongErr= 4.0533 
%
%coil 1 2 3 4  TrueErr= 4.5743  wrongErr=  4.6068
%
% concultion we fitting noise but maybe we can use voxel we not fitting to contorol it? but how?

%% let try a  nested CV fit

clist=[1 2 3];
   Fcoils{1}=[1 2 ]  ; Fcoils{2}=[1 3 ];Fcoils{3}=[2 3 ];
   Tcoils{1}=[3 ]  ; Tcoils{2}=[2 ];Tcoils{3}=[1 ];
   
[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinearLnO(par,M0SN(:,clist),OutPut.pBasis,nVoxels,Fcoils,Tcoils)...
         ,double(g0(:,clist)),[],[],options);
  
     G = OutPut.pBasis(:,:)*res1(:,:);

PDfit = zeros(nVoxels,1);
for ii=1:nVoxels
    PDfit(ii) = G(ii,:)' \ M0SN(ii,clist)';
end


BLSim = pdBiLinearFit_1(M0SN(:,clist), OutPut.pBasis(:,:), ...
    0, 1 ,0, PDfit(:), 1, Par(:,clist));
 M0Fit=BLSim.M0Fit

 
 



