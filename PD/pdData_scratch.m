%% Working out the multiple coil case.
%
% TODO
%   1.  Change to the Fourier basis
%   2.  Weight by signal to noise
%   3.  Try the statistics toolbox Robustfit method
%   4.  Check how Kendrick weights the polynomial basis and read on the web
%   about this
%   5.  Keep worrying about the polynomial order
%   6.  Check the Legendre polynomials, which are orthogonal, and see if
%   they help reduce the problem with the noise
%
%% If you are in the mrQ directory, run this to set the path
addpath(genpath(fullfile(mrqRootPath)));

%% Run the script for the pdPolyPhantomOrder
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 3;      % Second order is good for up to 5 samples
nSamples = 4;      % The box is -nSamples:nSamples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box location
oFlag = true;

printImages = false;
smoothkernel=[];
% This produces the key variables for comparing data and polynomial
% approximations. We will turn it into a function before long.
% Variables include M0S_v, pBasis, params, SZ
[OutPut] = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, oFlag);
% mrvNewGraphWin; imagesc(OutPut.pBasis);
% tmp = reshape(OutPut.pBasis,9,9,9,20);
% showMontage(tmp(:,:,:,1))

percentError = 100*OutPut.percentError;
fprintf('Polynomial approximation to the data (percent error): %0.4f\n',percentError)

%% Fit poly to Ratios
coilList = [1,2,3,4];
% coilList = [1,2,3,4,5,6,7];
% coilList = [7:-1:1];
% 
% coilList = 1:32;
% coilList = 1:2:32;
% coilList = 1:4:32;

%% First, try a pure simulation

% Fit the relative polynomial gain functions for the selected coils
% The 'est' parameter are the polynomial coefficients for each of the
% coils.
% This calculation uses the ratio of the coil data.
% The returned coil gains are specified up to an unknown scalar (they are
% relative).

M0pairs = [];
[polyRatio,M0pairs] = polyCreateRatio(OutPut.M0S_v(:,coilList), OutPut.pBasis);
estGainCoefficients = polySolveRatio(polyRatio,M0pairs);

% calculate PD fits error and param error in comper to the % the params derived from the phantom.
Res = polyRatioErr(estGainCoefficients, OutPut.params(:,coilList), OutPut.SZ(1:3), OutPut.pBasis);

% Show how well it does for pure simulation
estCoilGains  = OutPut.pBasis*Res.estGainParams';
trueCoilGains = OutPut.pBasis*Res.trueGainParams';
mrvNewGraphWin; plot(trueCoilGains(:),estCoilGains(:),'.')
identityLine(gca);
xlabel('Estimated gains'); ylabel('True gains');

%%  With real data, the fits go badly. 
% This is what we have to fix.
M0pairs = [];
[polyRatio, M0pairs] = polyCreateRatio(OutPut.M0_v(:,coilList), OutPut.pBasis);
estGainCoefficients  = polySolveRatio(polyRatio,M0pairs);
% s = svd(polyRatioMat);
% mrvNewGraphWin; plot(s)
% grid on

% calculate PD fits error and param error in comper to the % the params derived from the phantom.
Res = polyRatioErr(estGainCoefficients,OutPut.params(:,coilList), OutPut.SZ(1:3), OutPut.pBasis);

% Show how well it does for pure simulation
estCoilGains  = OutPut.pBasis*Res.estGainParams';
trueCoilGains = OutPut.pBasis*Res.trueGainParams';
mrvNewGraphWin; plot(trueCoilGains(:),estCoilGains(:),'.')
identityLine(gca);
xlabel('Estimated gains'); ylabel('True gains');

%%
showMontage(Res.PD)
%%
mrvNewGraphWin;
plot(M0pairs(:,1) ./ M0pairs(:,2), trueCoilGains - estCoilGains, '.')
xlabel('M0 ratios')
ylabel('Gain error');

%%
TG=reshape(trueCoilGains,[OutPut.SZ(1:3) 2]);
TR=TG(:,:,:,1)./TG(:,:,:,2); % the true ratio as estimate with poly with no noise
EG=reshape(estCoilGains,[OutPut.SZ(1:3) 2]);
ER=EG(:,:,:,1)./EG(:,:,:,2);% the ratio as  been fitted by our solotion
DR=(OutPut.M0(:,:,:,coilList(1))./OutPut.M0(:,:,:,coilList(2))); % the ratio of the raw data
showMontage(TR-DR);  title('no noise ratio- data ratio')
showMontage(TR-ER) ;  title('no noise ratio- estimated ratio')
showMontage(DR-ER);   title('data ratio - estimated ratio')
showMontage(TR-DR);  title('no noise ratio- data ratio')
showMontage(TR-ER) ;  title('no noise ratio- estimated ratio')
showMontage(1-ER./TR);   title('1- (no noise ratio / estimated ratio)')

%% The comparison in data land
% showMontage(M0(:,:,:,1))
% bar = reshape(trueCoilGains,11,11,11,3);
% showMontage(bar(:,:,:,1)); title('Correct')
% tmp = reshape(estCoilGains,11,11,11,3);
% showMontage(tmp(:,:,:,1)); title('Estimated')



%% let's visualized the fit in 3D space
 M0s=reshape(OutPut.M0S_v,OutPut.SZ);
 showMontage((M0s(:,:,:,1)-OutPut.M0(:,:,:,1))./OutPut.M0(:,:,:,1));colormap hot
 showMontage((M0s(:,:,:,2)-OutPut.M0(:,:,:,2))./OutPut.M0(:,:,:,2));colormap hot
 showMontage((M0s(:,:,:,3)-OutPut.M0(:,:,:,3))./OutPut.M0(:,:,:,3));colormap hot




%% To visualize the simulation versus the fits
%M0S = reshape(M0S_v,SZ);

% We also use M0S for some calculations below

% lets hold each time one dimation and look on the center  inplain
% plotRawandSimProfile(nCoils,M0,M0S,[1 1 1],10)
% plotRawandSimProfile(nCoils,M0,M0S,round(SZ(1:3)/2),11)
% plotRawandSimProfile(nCoils,M0,M0S,SZ(1:3),12)




%%

%% the real Data; 

% We also know that the data are similar to the simulations.  But, they are
% not exact.  When we send in the data, which are fit by the simulations to
% within about 1.2 percent, we don't get a good result.

[Res.est, Res.polyRatioMat] = polySolveRatio(M0_v(:,coilList),pBasis);

% calculate PD fits error and param error in comper to the % the params derived from the phantom.
[Res.CoilCoefErr, Res.PDerr,  Res.estMatrix,  Res.ParMatrix, ...
    Res.G, Res.M00, Res.PD, Res.PDspaceErr]= ...
    polyRatioErr(Res.est,params(:,coilList),pBasis);

Res=Res_D;




%%  make simulation with  noise and smooth in space

% do it on bigger size voulume and then crop (so the smooth will be homgenius in the relevant voulume)
noiseLevel=5;

 clear M0SNS M0SNS_v M0SN M0SN_v  st ed N rSizeDD nVoxelsDD pMatrixDD sDD
[pMatrixDD,sDD] = polyCreateMatrix(nSamples*10,2,3);
rSizeDD = length(sDD);
nVoxelsDD = rSizeDD^nDims;
M0SN_v = zeros(nVoxelsDD,nCoils);
for ii=1:nCoils
    M0SN_v(:,ii)= pMatrixDD*params(:,ii) +randn(nVoxelsDD,1)*noiseLevel;
end
M0SN=reshape(M0SN_v,[rSizeDD rSizeDD rSizeDD  SZ(4)]);
for i=1:nCoils
    A=smooth3(M0SN(:,:,:,i),'gaussian',[21 21 21]);
     %   A=smooth3(M0SNSd(:,:,:,i));
    M0SNS(:,:,:,i)=A;                 
end
N=(rSizeDD-rSize)/2;
st=N+1;
ed=rSizeDD-N;

M0SNS=M0SNS(st:ed,st:ed,st:ed,:);
M0SN=M0SN(st:ed,st:ed,st:ed,:);
showMontage(M0S-M0SNS);title('sim- noisey sim smooth')
showMontage(-M0S-M0SN) ;title('sim- noisey sim')

M0SN_v=reshape(M0SN,prod(SZ(1:3)),SZ(4));
M0SNS_v=reshape(M0SNS,prod(SZ(1:3)),SZ(4));
clear st ed N rSizeDD nVoxelsDD pMatrixDD sDD
%% Fits and plots
coilList = [5 6];


%simultions
[ Res_S ]   = fitRatioandPlotPD(coilList,M0S_v,M0S,pMatrix,params, 'Sim',0);

%simulted noise
SNR = 20*log10(mean(mean(M0S_v(:,coilList))) /noiseLevel)
[ Res_SN]   = fitRatioandPlotPD(coilList,M0SN_v,M0SN,pMatrix,params, 'SimNoise',0);
% Smooth simulted noise
[ Res_SNS]   = fitRatioandPlotPD(coilList,M0SNS_v,M0SNS,pMatrix,params, 'SimNoiseSmooth',1);


% Data
[ Res_data]   = fitRatioandPlotPD(coilList,M0_v,M0,pMatrix,params, 'data',0);




%%

 showMontage(100* ((Res_S.R- Res_SN.R)./ Res_S.R ) ,[],[],[],[],7 );colormap hot;title('100X(RS- RSN)/RS')
 showMontage( 100* ((Res_S.R- Res_SNS.R)./ Res_S.R   )  ,[],[],[],[],8 );colormap hot;title('100X(RS- RSNS)/RS')
 showMontage( 100* ((Res_S.R- Res_data.R)./ Res_S.R )  ,[],[],[],[],9 );colormap hot;title('RS- Rdata/RS')

%%
G1=reshape(Res_S.G,[SZ(1:3) 2]);
G2=reshape(Res_SN.G,[SZ(1:3) 2]);

         plotRawandSimProfile(2,G1,G2,[1 1 1],10,{'no noise' 'noise'})


%%
% i should study if the noise effect is depend on the R are there R that
% are less susptipble to noise the other like when it far from 1? 
%noise as function of R

% i should try and see of constrain like the same Z variation. or miror
% image in X dimantion that decrice the serarch space can make the sulotion
% more robust.

%i need for the first point:
%a function that 1 allow me to ingect ratio instade of Mo. i will use that to gradually and locally change the noise.
%in respect to R or space or M0. of gain params
%
%for the secound point i need a way to reduce the number of parameters or
%the dimantion in the matrix. for that i need to write the eqations. and
%see what fulling .
% i also need to explore the data and see what is the resons of simitry.
% are the in all the D parameters or becuse of a joinning of them.




%% % can i smoth R will this help?

R1=smooth3(Res_SimNoise.R,'gaussian');

showMontage(100* ((Res_Sim.R- R1)./ Res_Sim.R ))


Rparams= polyfitPhantomCoef(Res_SimNoise.R(:),pMatrix);
R2=pMatrix*Rparams;
R2=reshape(R2,SZ(1:3));

degree=4
[Poly,str] = constructpolynomialmatrix3d(SZ(1:3),find(ones(SZ(1:3))),degree); 
          
                [Rparams] = fit3dpolynomialmodel(Res_SimNoise.R,logical(Res_SimNoise.R) ,degree);
             R2=Poly*Rparams';
R2=reshape(R2,SZ(1:3));
 showMontage(100* ((Res_Sim.R- R2)./ Res_Sim.R ));colormap hot; caxis([-1.5 1.5])

 % answer: no no really this is not diffrent then smooth the data.that was
 % not helpful