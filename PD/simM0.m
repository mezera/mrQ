function [M0SN, M0S,SNR, PD]=simM0(CoilGain,PD,noiseLevel,PDtype,plotResults)
%            [M0SN, M0S,SNR]=simM0(OutPut.M0S_v(:,coilList),noiseLevel,PD)
% This function simulate MRI M0 data with guession noise with noise Level
% (noiseLevel), given a coils gains and PD. MOSN= CoilGain*PD+ Noise
%
%  Input
%  CoilGain       the coil gain: Columns containing the 3D Gain values in a vector  (nVoxels x nCoils)
%  PD                a PD values that will multipal the Coil Gain  (nVoxels x 1)
%  noiseLevel    the leven of the normal noise that will be added to the signal
% PDtype          to estimate a specail PD types :
%                           (1) 'dots'        - every tenth PD eqal 10 all other PD eqal 1
%                           (2) 'lowfreq'   - PD= sin(0:pi);
%                           (3) 'highfreq' - PD= sin(0:10*pi);
%  plotResults   when true (1)  plot M0S vs MOSN

%
% OutPut:
% M0S             the simulated M0  MOS= CoilGain*PD (nVoxels x nCoils)
% M0SN          the simulated M0  MOSN= CoilGain*PD+ Noise (nVoxels x nCoils)
% SNR             the simulated Signal to Noise ratio in desibals
% PD               the simulated PD
%   AM/BW Copyright VISTASOFT Team 2013

%% define the the problam dimation and parameters
nCoils=size(CoilGain,2);
nVoxels=size(CoilGain,1);

if notDefined('noiseLevel')
    noiseLevel = 0;
end
if notDefined('plotResults')
    plotResults = 0;
end

%% define PD
if notDefined('PD')
    PD = ones(nVoxels,1);
end

if notDefined('PDtype')
else
    PD = ones(nVoxels,1);
    switch PDtype
        case {'dots', '1'}  %pick point
            PD(1:ceil(nVoxels/10):nVoxels)=10;
        case {'highfreq', '2'}
            x=linspace(0,pi,nVoxels*10);
            PD=sin(x);
        case {'lowFreq', '3'}
            x=linspace(0,pi,nVoxels);
            PD=sin(x);
        otherwise
            error('PDtype %d not built',PDtype);
    end
end

%% simulate M0
M0S=zeros(nVoxels,nCoils);
M0SN=zeros(nVoxels,nCoils);
for ii=1:nCoils
    M0S(:,ii)=CoilGain(:,ii).*PD; %simulate M0 -->M0S= G*PD
    SNR(ii) = 20*log10(mean(M0S(:,ii))/noiseLevel);
    M0SN(:,ii)= M0S(:,ii) +randn(nVoxels,1)*noiseLevel;   % add noise --> M0SN= M0S +Noise
end

%% plot
fprintf(' the simulation SNR: %0.4f\n',SNR);

if plotResults==1
    x=length(M0S(:));
mrvNewGraphWin
plot(1:x,M0S(:),'.-k',1:x,M0SN(:),'.-r')
legend('M0S','M0SN')
xlabel('voxels')
xlabel('M0')
end

end