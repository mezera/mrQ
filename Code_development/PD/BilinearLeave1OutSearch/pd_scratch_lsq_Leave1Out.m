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

%%
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
