function [PD, R1]=mrQ_simulate_PD(PDtype,nVoxels)
% simulate differnt type of PD
%mrQ_simulate_PD(PDtype,nVoxels)
% inputs:
% PDtype type name
%nVoxels number of voxels
PD=[]; R1=[];
eSize = round(nVoxels^.333);
cPos = round(eSize/2);

PDtype = mrvParamFormat(PDtype);

switch PDtype
    case {'phantom'}
        PD = ones(nVoxels,1);
        R1=PD* 1.75;
        R1=R1./1000;
    case {'singlepoint'}
        PD = 0.5*ones(nVoxels,1);
        PD(sub2ind([eSize,eSize,eSize],cPos,cPos,cPos)) = 1;
        % showMontage(reshape(PD,eSize,eSize,eSize));
    case {'smallregion'}
        PD = 0.5*ones(nVoxels,1);
        if eSize<=3
            [X,Y,Z]=meshgrid(cPos,cPos-1:cPos+1,cPos);
        else
            [X,Y,Z]=meshgrid(cPos-1:cPos+1,cPos-1:cPos+1, cPos-1:cPos+1);
            
        end
        PD(sub2ind([eSize,eSize,eSize],X(:),Y(:),Z(:))) = 1;
    case {'linearslope'}
        PD = zeros(eSize,eSize,eSize);
        for ii=1:eSize
            PD(:,:,ii) = 0.5*ii/eSize;
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
    case {'tissue1', '4'} %tissue
        PD = zeros(eSize,eSize,eSize);
        mask = zeros(eSize,eSize,eSize);
        
        WM=0.7;
        GM=0.85;
        CSF=0.95;
        loc=randperm(nVoxels);
        jump=round(nVoxels./3);
        
        loc1=(1:jump);
        PD(loc(loc1))=WM+ 0.1*randn(length(loc1),1);
        mask(loc(loc1))=1;
        R1(loc(loc1)) = (2.5./PD(loc(loc1))) - 0.95;
        
        loc1=(jump+1:2*jump);
        PD(loc(loc1))=GM+ 0.1*randn(length(loc1),1);
        mask(loc(loc1))=2;
        R1(loc(loc1)) = (2../PD(loc(loc1))) - 0.9;
        
        loc1=(2*jump+1:nVoxels);
        PD(loc(loc1))=CSF+  0.1*randn(length(loc1),1);
        mask(loc(loc1))=3;
        R1(loc(loc1)) = (2.9./PD(loc(loc1))) - 0.99;
        
        PD = PD(:);
        R1=R1(:)./1000;
    case {'tissue2', '5'} % tissue with slope
        PD = zeros(eSize,eSize,eSize);
        mask = zeros(eSize,eSize,eSize);
        
        WM=0.7;
        GM=0.85;
        CSF=0.95;
        loc=randperm(nVoxels);
        jump=round(nVoxels./3);
        
        loc1=(1:jump);
        PD(loc(loc1))=WM;
        mask(loc(loc1))=1;
        R1(loc(loc1)) = (2.5./PD(loc(loc1))) - 0.95;
        
        loc1=(jump+1:2*jump);
        PD(loc(loc1))=GM;
        mask(loc(loc1))=2;
        R1(loc(loc1)) = (2../PD(loc(loc1))) - 0.9;
        
        loc1=(2*jump+1:nVoxels);
        PD(loc(loc1))=CSF;
        mask(loc(loc1))=3;
        R1(loc(loc1)) = (2.9./PD(loc(loc1))) - 0.99;
        
        
        for ii=1:eSize
            PD(:,:,ii) =PD(:,:,ii) +  0.1*ii/eSize;
        end
        PD = PD(:);
        R1=R1(:)./1000;
        
    otherwise
        error('PDtype %d not built',PDtype);
end

% unless define diferently simulate a linear relation of R1 and PD
if (notDefined('R1') || isempty(R1))
    R1 = (2.5./PD) - 2.26;
    R1=R1./1000;
end



end