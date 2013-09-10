%%% Figure 3
%
% Illustrate how PD can be mesure from noise M0 while estimate the coil gain using
%  bilinear sulotions and regularization by T1
%
% The problem we pace is the effect of noise on this solotions
%
% AM/BW Vistaosft Team, 2013

%%  Make sure mrQ is on the path
addpath(genpath(fullfile(mrqRootPath)));
rmpath('folderName') 
%% Generate example parameters for the coils from the phantom data

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box 
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations. We
% will turn it into a function before long. The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);

boxSize = repmat(phantomP.rSize,1,nDims);

%% simulte PD
[X,Y Z] = meshgrid(-nSamples:nSamples,-nSamples:nSamples, -nSamples:nSamples);
R  = sqrt(X.^2 + Y.^2 + + Z.^2);
PD = sin(R)./R;PD(isnan(PD))=1;
PD = abs(PD);
PD = sqrt(sqrt(sqrt(PD)));



%% simultae coil Gain (we are using the poylnomyal fits to the phantom data a typical coil function)
%select a set of coils

%we can sort coils by minmial correlation between them

%             c=nchoosek(1:16,k);
%              Cor=ones(4,size(c,1))*100;
%             for kk=1:size(c,1)
%                 A=(corrcoef(phantomP.M0_v(:,c(kk,:))));
%                 %    Cor(k,kk)=(sum(A(:))-k)/2;
%                 Cor(k,kk)=(sum(abs(A(:)))-k)/((length(A(:))-k)*0.5);
%                 
%                 Cloc(k,kk,1:k)=c(kk,:);
%             end
%         [v ind]=sort(Cor(:)); %sort the coils by minimum corralation
%         [xx yy]=ind2sub(size(Cor),ind(1)); % find the combination with minimal corralation
%        coils=[squeeze(Cloc(xx,yy,1:xx))']


coils = [1 3 5 9];
% get those coil poylnomyal coeficents 
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multipal the polynomyals coeficents
% and the polynomyal basis.
G = phantomP.pBasis*GainPolyPar;


%% simulte MRI SPGR  signal with and with out Noise
noiseLevel = 2;   % ?? Units???
% simultate the M0 and T1 fits of multi SPGR images. 
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,true);
% MR_Sim is a stracture with multipal fielled that inculde the simulation
% inputs MR sigunalinputsand the calculate of this signal after
% fitting the signal eqation.

%%

NL   = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,phantomP.pBasis);


%% T1 reg
% find the best cross validation
kFold=2;
lambda1= [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0];
R1basis(1:phantomP.nVoxels,1) = 1;  R1basis(:,2) =MR_Sim.R1Fit;
R1basis=double(R1basis);
[X_valdationErr,   gEstT, resnorm, FitT useX, kFold ]=pdX_valdationLoop_2( lambda1,kFold,MR_Sim.M0SN,phantomP.pBasis,R1basis,[],[],[]);

% find the lamda that best to cross validate (minimal RMSE error)
BestReg = find(X_valdationErr(2,:)==min(X_valdationErr(2,:)))% 

%        figure;  semilogy(lambda1,X_valdationErr(2,:),'*-');         xlabel('lambda');ylabel('CV error');X_valdationErr(2,:)./min(X_valdationErr(2,:))
% use the best CV lamda to fit the data
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ]=pdCoilSearch_T1reg( lambda1(BestReg),MR_Sim.M0SN,phantomP.pBasis,R1basis,gEstT(:,:,1,BestReg));

% figure;plot(PD(:)./mean(PD(:)), NL_T1reg.PD(:)./mean(NL_T1reg.PD(:)),'*')

% save the fitting outputs 

%%
% PD with T1 reg
PD_T1reg = reshape(NL_T1reg.PD, boxSize);
scale       = PD(1,1,1)/PD_T1reg(1,1,1);
PD_T1reg  = PD_T1reg.*scale;

% PD with out T1 reg
PD_Noreg = reshape(NL.PD, boxSize);
scale       = PD(1,1,1)/PD_Noreg(1,1,1);
PD_Noreg  = PD_Noreg.*scale;
return

%%  PD estimation with T1 reg vs without 

mrvNewGraphWin;

MM = minmax([PD_T1reg PD PD_Noreg]);
hold on
plot(PD_Noreg(:),PD(:),'o' ,'MarkerSize',10,'MarkerFaceColor','b')
plot(PD_T1reg(:),PD(:),'or','MarkerSize',10)

xlabel('Estimated PD'); ylabel('True PD');
identityLine(gca);
xlim([MM(1) MM(2)]);ylim([MM(1) MM(2)]);
axis image; axis square
legend('PD estimate with T1 reg','PD estimate without T1 reg','Location','NorthWest')



%% ridge reg
kFold=2;
lambda1= [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5  1e0 0.5 1e-1 0];
[X_valdationErr,   gEstT, resnorm, FitT useX, kFold ]=pdX_valdationRidgeLoop( lambda1,kFold,MR_Sim.M0SN,phantomP.pBasis);
%        figure;  semilogy(lambda1,X_valdationErr(2,:),'*-');         xlabel('lambda');ylabel('CV error');X_valdationErr(2,:)./min(X_valdationErr(2,:))

BestReg = find(X_valdationErr(2,:)==min(X_valdationErr(2,:)))% 


NL_Ridge = pdBiLinearFit_lsqRidgeSeach(MR_Sim.M0SN,phantomP.pBasis,gEstT(:,:,1,BestReg),[],lambda1(BestReg));

%%
