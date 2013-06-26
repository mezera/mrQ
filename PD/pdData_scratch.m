%% Working out the multiple coil case.
%
%
% TODO
%   1.  Realistic gains if possible -done X
%   2.  Extend to 3D  X
%   3.  Add a fourth coil or even make a function for N-coils  X
%   4.  Summarize the error distribution maybe in PD space? Coil space?  X
%   5.  Virtual coil analysis
%   6.  What else?
%   7.  Make the solution with eig and stuff a function  X
%   8.  Make the printing out and comparison a simple function  X
%  9.  Realistic noise
%
 addpath(genpath('/home/avivm/mrQ'));


%% Make sure mrQ is on your path
%addpath(genpath('/home/avivm/mrQ/PD'));
%To Set up parameters for N realistic coils
nCoils = 32;
nDims  = 3;
pOrder = 2;
nSamples=3;
%%
%the real data
% get M0 real sample 
%4D
[M0 SZ meanVal ]= phantomGetData(nSamples,3);
% resahave to 2D
M0_v=reshape(M0,prod(SZ(1:3)),SZ(4));

%% this is phantom data so we can simulte it by polyinomyals
[pMatrix,s] = polyCreateMatrix(nSamples,2,3);
rSize = length(s);
nVoxels = rSize^nDims;

% let fit
% fit the phantom poly coef assumint the phantoms PD eqal ones
for i=1:SZ(4)
params(:,i)= polyfitPhantomCoef(M0_v(:,i),pMatrix);
end


M0S_v = zeros(nVoxels,nCoils);
for ii=1:nCoils
    M0S_v(:,ii)= pMatrix*params(:,ii);
end
% 
%we can make the 4D simulation data
M0S=reshape(M0S_v,SZ);

%% to visuralized the simulation versase the fits
% lets hold each time one dimation and look on the center  inplain
         plotRawandSimProfile(nCoils,M0,M0S,[1 1 1],10)
         plotRawandSimProfile(nCoils,M0,M0S,round(SZ(1:3)/2),11)
         plotRawandSimProfile(nCoils,M0,M0S,SZ(1:3),12)


%%  make simuation with  noise 

noiseLevel = 5;
for i=1:nCoils
    M0SN_v(:,i) = M0S_v(:,i) + randn(size(M0S_v,1),1)*noiseLevel;
end
M0SN=reshape(M0SN_v,SZ);
% This is the signal to noise in decibels


%%  smooth the noise in coils space

for i=1:nCoils
    A=smooth3(M0SN(:,:,:,i),'gaussian');
    M0SNs(:,:,:,i)=A;
end
M0SNs_v=reshape(M0SNs,[prod(SZ(1:3)) SZ(4)]);
% get the gain estimates for the noisy data


%% lets add the noise ans smoth on a double size box and then crop 
%(so the smooth will be homgenius in the ROI)
noiseLevel=5;


[pMatrixDD,sDD] = polyCreateMatrix(nSamples*2,2,3);
rSizeDD = length(sDD);
nVoxelsDD = rSizeDD^nDims;
M0SNSd_v = zeros(nVoxelsDD,nCoils);
for ii=1:nCoils
    M0SNSd_v(:,ii)= pMatrixDD*params(:,ii) ++ randn(nVoxelsDD,1)*noiseLevel;
end
M0SNSd=reshape(M0SNSd_v,[(nSamples*4+1) (nSamples*4+1) (nSamples*4+1)  SZ(4)]);
for i=1:nCoils
%    A=smooth3(M0SNSd(:,:,:,i),'gaussian',[9 9 9]);
        A=smooth3(M0SNSd(:,:,:,i));

    M0SNSd(:,:,:,i)=A;                 
end
st=nSamples*2+1-nSamples;
ed=nSamples*2+1+nSamples;

M0SNSd=M0SNSd(st:ed,st:ed,st:ed,:);

M0SNSd_v=reshape(M0SNSd,[prod(SZ(1:3)) SZ(4)]);
clear  st ed A pMatrixDD rSizeDD sDD nVoxelsDD


%% Fits and plots
coilList = [1 2];


%simultions
[ Res_Sim ]   = fitRatioandPlotPD(coilList,M0S_v,M0S,pMatrix,params, 'Sim',0);

%simulted noise
SNR = 20*log10(mean(mean(M0S_v(:,coilList))) /noiseLevel)
[ Res_SimNoise]   = fitRatioandPlotPD(coilList,M0SN_v,M0SN,pMatrix,params, 'SimNoise',0);
% Smooth simulted noise
[ Res_SimNoiseSH1]   = fitRatioandPlotPD(coilList,M0SNs_v,M0SNs,pMatrix,params, 'SimnoiseSH1',0);
[ Res_SimNoiseSH2]   = fitRatioandPlotPD(coilList,M0SNSd_v,M0SNSd,pMatrix,params, 'SimnoiseSH2',0);


% Data
[ Res_data]   = fitRatioandPlotPD(coilList,M0_v,M0,pMatrix,params, 'data',0);




%%

 showMontage(100* ((Res_Sim.R- Res_SimNoise.R)./ Res_Sim.R ) ,[],[],[],[],7 );colormap hot;title('Rsim- RsimNoise/Rsim')
 showMontage( 100* ((Res_Sim.R- Res_SimNoiseSH1.R)./ Res_Sim.R  ) ,[],[],[],[],8  );colormap hot;title('Rsim- RsimNoiseSH1/Rsim')
 showMontage( 100* ((Res_Sim.R- Res_SimNoiseSH2.R)./ Res_Sim.R   )  ,[],[],[],[],9 );colormap hot;title('Rsim- RsimNoiseSH2/Rsim')
 showMontage( 100* ((Res_Sim.R- Res_data.R)./ Res_Sim.R )  ,[],[],[],[],10 );colormap hot;title('Rsim- Rdata/Rsim')


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