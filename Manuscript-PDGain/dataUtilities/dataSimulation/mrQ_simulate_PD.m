function [PD, R1]=mrQ_simulate_PD(PDtype,nVoxels)
% Simulate different volume patterns of PD and R1
%
%    mrQ_simulate_PD(PDtype,nVoxels)
% 
% Inputs:  
%  PDtype  - Type of PD structure in the simulation 
%    0 - Uniform
%    1 - Dot arrays
%    2 - High freq
%    3 - Low frequency
%    4,5,6 - Tissue types
%    7 - Single point
%    8 - Small Region
%    9 - Linear slope
%
%  nVoxels - Total number of voxels in the volume.  It is adjusted to be a
%  cube, though
%
%
% AM Vistasoft Team, 2013

%% Initialize the parameters

% The edge size of the box
eSize = round(nVoxels^(1/3));

% Center position
cPos = round(eSize/2);

% This is the parameter we have been using to specify the box size
nSamples=(eSize-1)/2;

%% Here is the flag.  Let's go
PDtype = mrvParamFormat(PDtype);
switch PDtype
    case {'phantom' '0'}
        % Uniform
        PD = ones(nVoxels,1);
        R1 = PD* 1.75;
        R1 = R1./1000;
        PD = PD(:);
    case {'dots', '1'}  
        % BROKEN
        PD = ones(nVoxels,1);
        PD(1:ceil(nVoxels/10):nVoxels) = 10;
    case {'highfreq', '2'}
        x = linspace(0,pi*10,nVoxels);
        PD = sin(x);
    case {'lowfreq', '3'}
        x = linspace(0,pi,nVoxels);
        PD = sin(x);
        disp('Low frequency')
    case {'tissue1', '4'} 
        % Tissue properties simulated for gray, white and CSF
        PD = zeros(eSize,eSize,eSize);
        % mask = zeros(eSize,eSize,eSize);
        
        WM = 0.7;    % Water fraction in the different tissue types
        GM = 0.85;
        CSF= 0.95;
        
        loc = randperm(nVoxels);
        jump= round(nVoxels./3);
        
        % Start with R1 units in 1/sec
        loc1=(1:jump);
        PD(loc(loc1))  = WM + 0.1*randn(length(loc1),1);
        % mask(loc(loc1))= 1;
        R1(loc(loc1))  = (2.5./PD(loc(loc1))) - 0.95;
        
        loc1=(jump+1:2*jump);
        PD(loc(loc1))  = GM+ 0.1*randn(length(loc1),1);
        % mask(loc(loc1))= 2;
        R1(loc(loc1))  = (2../PD(loc(loc1))) - 0.9;
        
        loc1=(2*jump+1:nVoxels);
        PD(loc(loc1))  = CSF+  0.1*randn(length(loc1),1);
        % mask(loc(loc1))= 3;
        R1(loc(loc1))  = (2.9./PD(loc(loc1))) - 0.99;
        
        PD = PD(:);
        R1 = R1(:)./1000;    % Convert units to milliseconds
    case {'tissue2', '5'} 
        % Tissue with slight PD slope and deviation in T1-PD relationship
        PD = zeros(eSize,eSize,eSize);
        % mask = zeros(eSize,eSize,eSize);
        
        WM=0.7;  GM = 0.85;  CSF=0.95;
        loc  = randperm(nVoxels);
        jump = round(nVoxels./3);
        
        loc1 = (1:jump);
        PD(loc(loc1))   = WM;
        % mask(loc(loc1)) = 1;
        R1(loc(loc1))   = (2.5./PD(loc(loc1))) - 0.95;
        
        loc1=(jump+1:2*jump);
        PD(loc(loc1))   = GM;
        % mask(loc(loc1)) = 2;
        R1(loc(loc1))   = (2../PD(loc(loc1))) - 0.9;
        
        loc1=(2*jump+1:nVoxels);
        PD(loc(loc1))   = CSF;
        % mask(loc(loc1)) = 3;
        R1(loc(loc1))   = (2.9./PD(loc(loc1))) - 0.99;
        
        % Put a little spatial structure in the PD and break the PD-T1
        % relationship a bit
        for ii=1:eSize
            PD(:,:,ii) =PD(:,:,ii) +  0.1*ii/eSize;
        end
        PD  = PD(:);
        R1  = R1(:)./1000;
        
    case {'tissue3', '6'}
        % Spatial pattern in the PD
        [X,Y, Z] = meshgrid(-nSamples:nSamples,-nSamples:nSamples, -nSamples:nSamples);
        R  = sqrt(X.^2 + Y.^2 + Z.^2);
        
        % R is the distance from the center.  We make a rectified sinusoid
        % from the center to the edge.  We set all the NaN values to 1.  We
        % then take the sixth root to squeeze the dynamic range to be
        % reasonable.
        PD            = sin(R)./R;
        PD(isnan(PD) )= 1;
        PD = abs(PD);
        PD = PD .^ (1/6);
        
    case {'singlepoint','7'}
        % PD = 0.5 everywhere except one point where it is 1
        PD = 0.5*ones(nVoxels,1);
        PD(sub2ind([eSize,eSize,eSize],cPos,cPos,cPos)) = 1;
        % showMontage(reshape(PD,eSize,eSize,eSize));
    case {'smallregion','8'}
        % PD = 0.5 except in a small region
        PD = 0.5*ones(nVoxels,1);
        if eSize<=3
            [X,Y,Z]=meshgrid(cPos,cPos-1:cPos+1,cPos);
        else
            [X,Y,Z]=meshgrid(cPos-1:cPos+1,cPos-1:cPos+1, cPos-1:cPos+1);
            
        end
        PD(sub2ind([eSize,eSize,eSize],X(:),Y(:),Z(:))) = 1;
        
    case {'linearslope','9'}
        % Linear slope of the PD
        PD = zeros(eSize,eSize,eSize);
        for ii=1:eSize
            PD(:,:,ii) = 0.5*ii/eSize;
        end   
        
    otherwise
        error('PDtype %d not built',PDtype);
end

% Convert to a volume
PD = reshape(PD,eSize,eSize,eSize);

% Unless defined above, simulate a linear relation of R1 and PD
if (notDefined('R1') || isempty(R1))
    R1 = (2.5./PD) - 2.26;
    R1 = R1./1000;
end


end