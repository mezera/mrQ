function [Boxes, PositiveBoxs, UnCorBoxs, UnSTDBoxs]=mrQ_CalBoxPD_step1a(opt,BoxesToUse,CoilGains)

% [Boxes, PositiveBoxs, UnCorBoxs, UnSTDBoxs]=mrQ_CalBoxPD_step1a(opt, ...
%                                                     BoxesToUse,CoilGains)
% This is step 1 of 6 for building the WF map.
%     Step 1: Get the x, y, z, and PD values of the boxes. 
%
% ~INPUTS~ 
%          opt:
%   BoxesToUse:
%    CoilGains:
%
% ~OUTPUTS~
%        Boxes:
% PositiveBoxs:
%    UnCorBoxs:
%    UnSTDBoxs:
%
% See also: mrQ_buildPD_ver2
%           Step_0: none
%           Step_2: mrQ_ScaleBoxes_step2
%           Step_3: mrQ_BoxJoinBox
%           Step_4: mrQ_smoothGain_step4b
%           Step_5: mrQ_PD2WF_step5
%
% AM Vistalab team 2013


%% I. Book keeping

% bookkeeping of the box fits that went wrong
PositiveBoxs=zeros(length(opt.wh),1);
UnCorBoxs=zeros(length(opt.wh),1);
UnSTDBoxs=zeros(length(opt.wh),1);

%% II. Get M0, T1 and Brain Mask 
% Get the M0 and T1 information
% multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
SZ=size(M0);
%T1
%T1=readFileNifti(opt.T1file);
%T1=T1.data;

%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

smoothkernel=opt.smoothkernel;
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);
maxCoil=opt.maxCoil;
minCoil=opt.minCoil;
useCoil=opt.useCoil;

if length(SZ)==4;        
    if length(useCoil)>SZ(4); 
        useCoil=useCoil(1:SZ(4));
    end; 
end

nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

Ncoils=length(useCoil);

for ii=BoxesToUse
    clear M01  t1  BM1  SZ M0_v R1basis PD
    
    [fb(1,1,1), fb(1,1,2), fb(1,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, ~, BM1, SZ, skip, ~,~, XX, YY, ZZ ]= mrQ_GetM0_boxData(opt,[],M0,BM,fb(1,1,:),smoothkernel);
    if length(SZ)==4;
    M0_v = reshape(M01, prod(SZ(1:3)), SZ(4));
    else
        M0_v=M01(:);
    end
        
    %make a R1 regularization matrix
  %  R1basis(1:nVoxels,1) = 1; R1=1./(t1(:)); R1basis(:,2) =R1;
  %  R1basis=double(R1basis);

%% III. G, PD, coil gain and coil PD
    % 1. Get G
    g=CoilGains(ii).g;
    G = pBasis*g;
    
    % 2. Solve PD
    Clist=CoilGains(ii).Clist;
    PD = zeros(nVoxels,1);
    for jj=1:nVoxels
        PD(jj) = G(jj,:)' \ M0_v(jj,Clist)';
    end
    
% If there are zeroes, we will get NaN. 
% If it less then zero, it is just wrong
    mask=PD>0; 
    
% Check that the between-coil error is minimal
    PDC=M0_v(:,Clist)./G;
    EstimateErr=std(PDC./repmat(mean(PDC,2),1,size(PDC,2)),[],2);
    mask=mask & EstimateErr<0.08;
    
           PDSTD=std(PD(mask));
           
% Check that there are no outlier voxels
    mask=mask & PD<PD+3*PDSTD & PD>PD-3*PDSTD ;
    
  
    % 3. Solve for all coils Gain
    G  = zeros(nVoxels,Ncoils);
    g0 = zeros(nPolyCoef,Ncoils);
    for jj=1:Ncoils
        G(mask,jj)  = M0_v(mask,jj) ./ PD(mask);         % Raw estimate
        g0(:,jj) = pBasis(mask,:) \ G(mask,jj);  % Polynomial approximation
    end
    G = pBasis*g0;
    
    % 4. Solve all coils PD
    PD = zeros(nVoxels,1);
    V=1:nVoxels & mask';
    for jj=find(V)
        PD(jj) = G(jj,:)' \ M0_v(jj,useCoil)';
    end
    
    Boxes(ii).PD=PD;
    Boxes(ii).XX=XX;
    Boxes(ii).YY=YY;
    Boxes(ii).ZZ=ZZ;
 
%% IV. Various Checks    
% Check for correlation role
    if length(SZ)==4
    Mc=corrcoef(M0_v(:,Clist));Gc=corrcoef(G(:,Clist));
    Bad= find(Mc<Gc);
    else
        Bad=[];
    end
    
% Check for STD error
if length(SZ)==4
    W=mean(M0_v(mask,useCoil),2)./Ncoils;
    Dat=sum(M0_v(mask,useCoil).*repmat(W,1,Ncoils),2);
    Dat=(Dat./mean(Dat(:))).*mean(PD(mask));
    DatSTD=std(Dat);
    PDSTD=std(PD(mask));
    
    if DatSTD*1.1<PDSTD
        UnSTDBoxs(ii)=1;
    end
else
    UnSTDBoxs(ii)=0;
end
    
% Check for bad solution negative PD. 
    if     length(find(mask))/length(mask)<0.5 
%        If there are 50% negative values, it is clearly a wrong solution
         Boxes(ii).NegativeBad=1;
       
         
    else
         Boxes(ii).NegativeBad=0;
         PositiveBoxs(ii)=1;
         
         mask1=zeros(size(BM));
         mask1(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2))=1;
         wh=find(mask1);
         Boxes(ii).loc=wh(mask);
          Boxes(ii).PD=PD(mask);
          
    end
    
% Check if the M0 vector is more correlated than the G. 
    % They must be or the solution is wrong.
    if isempty(Bad)
        Boxes(ii).Corgood=1;
        UnCorBoxs(ii)=1;
    else
        Boxes(ii).Corgood=0;
    end
    Boxes(ii).CorBad=Bad;
    
    
end
