function [M0SN, M0S, SNR, PD]=simM0(CoilGain, PD, noiseLevel, plotFlag)
% Simulate M0 data with gaussian noise
%
% [M0SN, M0S,SNR] =simM0(OutPut.M0S_v(:,coilList), PD, noiseLevel, plotFlag)
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
%  plotFlag when true (1) plot M0S vs MOSN
%
% OutPut:
%  M0S     Noise free simulated M0  MOS= CoilGain*PD (nVoxels x nCoils)
%  M0SN    The simulated M0  MOSN= CoilGain*PD+ Noise (nVoxels x nCoils)
%  SNR     The simulated Signal to Noise ratio in desibals
%  PD      The simulated PD
%
% AM/BW Copyright VISTASOFT Team 2013

%% Define the the problem dimension and parameters

nCoils  = size(CoilGain,2);
nVoxels = size(CoilGain,1);

if notDefined('noiseLevel'),  noiseLevel = 0; end
if notDefined('plotFlag'), plotFlag = 0; end

%% define PD
if notDefined('PD'), PD = ones(nVoxels,1); end
if ischar(PD)
    PD = mrvParamFormat(PD);
    eSize = round(nVoxels^.333);
    cPos = round(eSize/2);
    switch PD
        case {'singlepoint'}
            PD = 0.5*ones(nVoxels,1);
            PD(sub2ind([eSize,eSize,eSize],cPos,cPos,cPos)) = 1;
            % showMontage(reshape(PD,eSize,eSize,eSize));
        case {'smallregion'}
            PD = 0.5*ones(nVoxels,1);
            PD(sub2ind([eSize,eSize,eSize],X(:),Y(:),Z(:))) = 1;
        case {'linearslope'}
            PD = zeros(eSize,eSize,eSize);
            for ii=1:eSize
                PD(:,:,ii) = ii/eSize;
            end
            PD = PD(:);
        case {'dots', '1'}  %pick point
            PD = ones(nVoxels,1);
            PD(1:ceil(nVoxels/10):nVoxels) = 10;
        case {'highfreq', '2'}
            x = linspace(0,pi,nVoxels*10);
            PD = sin(x);
        case {'lowFreq', '3'}
            x = linspace(0,pi,nVoxels);
            PD = sin(x);
        otherwise
            error('PDtype %d not built',PDtype);
    end
end

%% simulate M0
M0S  = zeros(nVoxels,nCoils);   % Noise free
M0SN = zeros(nVoxels,nCoils);   % M0 plus gaussian noise
SNR  = zeros(nCoils,1);         % Save out the SNR
%
for ii=1:nCoils
    M0S(:,ii)=CoilGain(:,ii).*PD;     %simulate M0 -->M0S= G*PD
    SNR(ii) = 20*log10(mean(M0S(:,ii))/noiseLevel);  % Coil SNR
    % Simulate M0 plus noise: M0SN= M0S +Noise
    M0SN(:,ii)= M0S(:,ii) + randn(nVoxels,1)*noiseLevel;
end

%% plot
fprintf('The simulation SNR: %0.4f\n',SNR);

if plotFlag 
    mrvNewGraphWin

    x = length(M0S(:));
    
    plot(1:x,M0S(:),'.-k',1:x,M0SN(:),'.-r')
    legend('M0S','M0SN')
    xlabel('voxels')
    ylabel('M0')
end

end