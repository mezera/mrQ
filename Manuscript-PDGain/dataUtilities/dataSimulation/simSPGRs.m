function [OutPutSim]=simSPGRs(CoilGain, PD, B1,R1,flipAngles,tr, noiseLevel, plotFlag)
% Simulate SPGR signal and the M0 T1 fit with Gaussian noise
%
% [OutPutSim] = simSPGRs(CoilGain, PD, B1,R1,flipAngles,tr, noiseLevel, plotFlag)
%
% The noise model is
%
%      MOsimulated = (CoilGain * PD)  + Gaussian-Noise
%
% Input
%  CoilGain       Columns containing the 3D Gain values in a vector (nVoxels x nCoils)
%  PD             PD values that will multiply the Coil Gain  (nVoxels x 1)
%  noiseLevel     the level of the noise
%  PDtype         Special PD patterns:
%                    (1) 'dots'     - every tenth PD eqal 10 all other PD eqal 1
%                    (2) 'lowfreq'  - PD= sin(0:pi);
%                    (3) 'highfreq' - PD= sin(0:10*pi);
% plotFlag    When true (1) plot M0S vs MOSN
% B1          The excite in-homogenity at a number that multipal the
%              flip angel in each location (nVoxels x 1). (default no
%              in-homogenity (B1=1))
% R1          The R1 value is msec-1 in each location(nVoxels x 1). If not
%               defined the defult will be a linear function of PD.
% flipAngles  The flip angles in degrees. for more then one flipAngles
%               (flipAngles X 1). default flipangle=[4 10 20 30];
% tr          Time to repeat is msec defult (tr=20)
% te          Not yet implemented and assumed to be not important with
%                very short te
% OutPut:
%  M0S     Noise free simulated M0  MOS= CoilGain*PD (nVoxels x nCoils)
%  M0SN    The simulated M0  MOSN= CoilGain*PD+ Noise (nVoxels x nCoils)
%  SNR     The simulated Signal to Noise ratio in desibals
%  PD      The simulated PD
%  mask    Tissue mask in case of tissue simulation 1=WM GM=2 CSF=3;
%
% AM/BW Copyright VISTASOFT Team 2013

%% Define the the problem dimension and parameters
mask=[];
nCoils  = size(CoilGain,2);
nVoxels = size(CoilGain,1);

if notDefined('noiseLevel'),  noiseLevel = 0; end
if notDefined('plotFlag'), plotFlag = 0; end
if notDefined('flipAngles'), flipAngles=[4 10 20 30];end
if notDefined('B1'), B1( 1:nVoxels,1)= 1; end
if notDefined('tr'), tr= 20; end

%% define PD
if notDefined('PD'), PD = ones(nVoxels,1); end
if ischar(PD)
    PD = mrvParamFormat(PD);
    [PD, R11]=mrQ_simulate_PD(PD,nVoxels) 
    if (notDefined('R1') || isempty(R1))
        R1=R11;
    end
end
%when R1 is not defined or simulated we will use the literature linear
%relations with PD
if (notDefined('R1') || isempty(R1))
    R1 = (2.5./PD) - 2.26;
    R1=R1./1000;
end



%% simulate M0
M0S  = zeros(nVoxels,nCoils);   % M0 Noise free
Sig = zeros(nVoxels,length(flipAngles),nCoils);   % Sig Noise free
SigN = zeros(nVoxels,length(flipAngles),nCoils);   % Sig plus gaussian noise
SNR  = zeros(nCoils,length(flipAngles));         % Save out the SNR

%
%
% the SPGR eqation
for jj=1:length(flipAngles)
    for ii=1:nCoils
        %the affactive flip Angles is also depend on B1
        fa=flipAngles(jj).*B1;
        fa = fa./180.*pi;
        
        
        M0S(:,ii)=CoilGain(:,ii).*PD;     %simulate M0 -->M0S= G*PD
        Sig(:,jj,ii) =M0S(:,ii).*(1-exp(-tr.*R1)).*sin(fa)./(1-exp(-tr.*R1).*cos(fa));
        
        
        SNR(ii,jj) = 20*log10(mean(Sig(:,jj,ii))/noiseLevel);  % Coil SNR
        % Simulate M0 plus noise: M0SN= M0S +Noise
        SigN(:,jj,ii)=  Sig(:,jj,ii) + randn(nVoxels,1)*noiseLevel;
    end
end

%% Fit The R1 and M0
options = optimset('Algorithm', 'levenberg-marquardt','Display', 'off','Tolx',1e-12);
%get the sum of sqr for T1 and M0 fit
SIG=sqrt(sum(SigN.^2,3));
% initilize the fit by the linear aproxsimation of the signal eqation sulotion
[t1t,M0t]     = relaxFitT1(SIG,flipAngles,tr,B1);
% fit each voxel with lsq non linear fit on the full single eqation
for ii= 1:nVoxels
    x0(1)=M0t(ii);
    x0(2)=t1t(ii)*1000;
    [res(:,ii), resnorm(ii)] = lsqnonlin(@(par) errT1PD(par,flipAngles,tr,SIG(ii,:),1,1,1,[]),x0,[],[],options);
    
end

M0Fit  = res(1,:);
R1Fit  = 1./res(2,:);

%%  calcultate the fitted M0 per coil given the T1 fit
for jj=1:length(flipAngles)
    for ii=1:nCoils
        fa=flipAngles(jj).*B1;
        fa = fa./180.*pi;
        M0SNf(:,jj,ii)= (SigN(:,jj,ii).*(1-exp(-tr.*R1Fit(:)).*cos(fa)))./((1-exp(-tr.*R1Fit(:))).*sin(fa));
        
    end
end
M0SN=squeeze( mean(M0SNf,2));



%% OutPutSim
OutPutSim.M0SN=M0SN;
OutPutSim.M0S=M0S;
OutPutSim.SNR=SNR;
OutPutSim.PD=PD;
OutPutSim.mask=mask;
OutPutSim.R1=R1;
OutPutSim.R1Fit=R1Fit;
OutPutSim.SigN=SigN;
OutPutSim.Sig=Sig;
OutPutSim.tr=tr;
OutPutSim.flipAngles=flipAngles;


%% plot
if plotFlag
fprintf('The simulation SNR: %0.4f\n',SNR);
end
if plotFlag
    mrvNewGraphWin
    subplot(1,2,1)
    x = length(M0S(:));
    
    plot(1:x,M0S(:),'.-k',1:x,M0SN(:),'.r')
    legend('M0S','M0SN Fit')
    xlabel('voxels')
    ylabel('M0')
    
    subplot(1,2,2)
    x = length(R1(:));
    
    plot(1:x,R1(:),'.-k',1:x,R1Fit(:),'.r')
    legend('R1','R1 Fit')
    xlabel('voxels')
    ylabel('R1')
    
end

end