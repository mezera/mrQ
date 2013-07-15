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
Par=OutPut.params(:,[1 :8]);
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

%   mean of squr
 PDsosq = sqrt(sum(M0SN.^2,2));
 PDinit=PDsosq;
 PDinit=PDinit(:);

%   random
% PDinit = rand(size(PDsim(:)));
% PDinit=PDinit(:);

%   segmentaion
% PDinit=nan(size(mask));
%  PDinit(find(mask==1))=1;
% PDinit=PDinit(:);
%
%    start with  reage regration fit
% PDinit = rand(size(PDsim(:)));
% PDinit=PDinit(:);
% D = diag(OutPut.W); %D(1,1)=0.01;
% % diag(D)
% % Constant terme
% D(1,1) = 0;
% % Linear terms -
% lWeight = .01;
% sWeight = 1;
% cWeight = 1;
% D(2,2) = lWeight;  D(4,4) = lWeight; D(7,7) = lWeight;
% % Quadratic terms
% D(3,3) = sWeight; D(5,5)= sWeight; D(6,6) = sWeight;
% % Cross products?
% D(8,8) = cWeight; D(9,9) = cWeight; D(10,10) = cWeight;
% BLSim = pdBiLinearFit_1(M0SN, OutPut.pBasis, ...
%     1, 100, 0, PDinit(:), 1, Par,D);
%  PDinit = BLSim.PD(:);



%   true sulotiop
%PDinit=PDsim(:);



% get inital guess
% G = zeros(nVoxels,nCoilsS);
% g0 = zeros(nPolyCoef,nCoilsS);
% we can be spesipic with what we start the rest will be zeros.
mask1=~isnan(PDinit);
for ii=1:nCoilsS
    G(mask1,ii)  = M0SN(mask1,ii) ./ PDinit(mask1);         % Raw estimate
    g0(:,ii) =OutPut. pBasis(mask1,:) \ G(mask1,ii);  % Polynomial approximation
end


%% let try a  nested CV fit
options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',inf,'TolFun', 1e-8,'TolX', 1e-10);

clist=[1 2 3 ];
Fcoils{1}=[1 2 ]  ; Fcoils{2}=[1 3 ];Fcoils{3}=[2 3 ];
Tcoils{1}=[3 ]  ; Tcoils{2}=[2 ];Tcoils{3}=[1 ];



% clist=[1 2 3 4 5 6];
% Fcoils{1}=[1 2 ]  ; Fcoils{2}=[1 3 ];Fcoils{3}=[2 3 ];
% Tcoils{1}=[4 ]  ; Tcoils{2}=[5 ];Tcoils{3}=[6 ];

[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinearLnO(par,M0SN(:,clist),OutPut.pBasis,nVoxels,Fcoils,Tcoils)...
    ,double(g0(:,clist)),[],[],options);

G = OutPut.pBasis(:,:)*res1(:,:);

PDfit = zeros(nVoxels,1);
for ii=1:nVoxels
    PDfit(ii) = G(ii,:)' \ M0SN(ii,clist)';
end


BLSim = pdBiLinearFit_1(M0SN(:,clist), OutPut.pBasis(:,:), ...
    0, 1 ,0, PDinit(:), 1, Par(:,clist));
M0Fit=BLSim.M0Fit


RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
err= (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:)) );
err=reshape(err,OutPut.SZ(1:3));
showMontage(err);title(['RMSE : ' num2str(RMSE)]);
%% i got the feeling that my result "forget" to keep the PD constant i need to impose it
for jj=1:length(clist)
        PDj(:,jj) = M0SN(:,clist(jj))./G(:,clist(jj)) ;
end
figure;plot(std(PDj,[],2)./median(PDj,2))
    
%% stocastic sets



options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',5,'TolFun', 1e-6,'TolX', 1e-10);
PFfit0=zeros(nVoxels,1);
PDchange=[];
Fcoils{1}=[1 2]  ; Fcoils{2}=[1 3];Fcoils{3}=[2 3 ];

Tcoils{1}=3   ; Tcoils{2}=2; Tcoils{3}=1;

ii=1
while ii<100
    
    CC=randperm(nCoilsS);
    clist=CC((1:3))
    
    if ii==1
        g1=g0(:,clist);
    else
        G = zeros(nVoxels,3);
        g1 = zeros(nPolyCoef,3);
        %
        for jj=1:3
            G(:,jj)  = M0SN(:,clist(jj)) ./ PDfit(:);         % Raw estimate
            g1(:,jj) =OutPut. pBasis \ G(:,jj);  % Polynomial approximation
        end
        
    end
    
    [res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinearLnO(par,M0SN(:,clist),OutPut.pBasis,nVoxels,Fcoils,Tcoils)...
        ,double(g1),[],[],options);
    
    G = OutPut.pBasis*res1;
    PDfit = zeros(nVoxels,1);
    for jj=1:nVoxels
        PDfit(jj) = G(jj,:)' \ M0SN(jj,clist)';
    end
    PDchange(ii)=std(PDfit-PFfit0);
    PFfit0=PDfit;
    ii=ii+1
    
end


%%

BLSim = pdBiLinearFit_1(M0SN(:,clist), OutPut.pBasis(:,:), ...
    0, 1 ,0, PDfit(:), 1, Par(:,clist));
M0Fit=BLSim.M0Fit


RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
err= (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:)) );
err=reshape(err,OutPut.SZ(1:3));
showMontage(err);title(['RMSE : ' num2str(RMSE)]);

