function Boxes=mrQ_CalBoxPD_step1(opt,BoxesToUse,CoilGains)
% Boxes=mrQ_CalBoxPD_step1(opt,BoxesToUse,CoilGains)
%
%
% AM Vistalab team 2013

% Get the M0 and T1 information

%multi coil M0
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
        if length(useCoil)>SZ(4); useCoil=useCoil(1:SZ(4));end

nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

Ncoils=length(useCoil);

for ii=BoxesToUse
    clear M01  t1  BM1  SZ M0_v R1basis PD
    
    [fb(1,1,1), fb(1,1,2), fb(1,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, ~, BM1, SZ, skip, ~,~, XX, YY, ZZ ]= mrQ_GetM0_boxData(opt,[],M0,BM,fb(1,1,:),smoothkernel);
    
    M0_v = reshape(M01, prod(SZ(1:3)), SZ(4));
    %make a R1 regularazation matrix
  %  R1basis(1:nVoxels,1) = 1; R1=1./(t1(:)); R1basis(:,2) =R1;
  %  R1basis=double(R1basis);
    
    % 1. Get G
    g=CoilGains(ii).g;
    G = pBasis*g;
    
    % 2.  solve PD
    Clist=CoilGains(ii).Clist;
    PD = zeros(nVoxels,1);
    for jj=1:nVoxels
        PD(jj) = G(jj,:)' \ M0_v(jj,Clist)';
    end
    %if there are zeors we will get NAN. 
    %if it less then Zerow it is just wrong
    mask=PD>0; 
    
    % 3. solve all coils Gain
    G  = zeros(nVoxels,Ncoils);
    g0 = zeros(nPolyCoef,Ncoils);
    for jj=1:Ncoils
        G(mask,jj)  = M0_v(mask,jj) ./ PD(mask);         % Raw estimate
        g0(:,jj) = pBasis(mask,:) \ G(mask,jj);  % Polynomial approximation
    end
    G = pBasis*g0;
    % 4 solve all coil PD
    PD = zeros(nVoxels,1);
    for jj=1:nVoxels
        PD(jj) = G(jj,:)' \ M0_v(jj,useCoil)';
    end
    
    Boxes(ii).PD=PD;
    Boxes(ii).XX=XX;
    Boxes(ii).YY=YY;
    Boxes(ii).ZZ=ZZ;
    Mc=corrcoef(M0_v(:,useCoil));Gc=corrcoef(G);
    Bad= find(Mc<Gc);
    if isempty(Bad)
        Boxes(ii).Corgood=1;
    else
        Boxes(ii).Corgood=0;
    end
    Boxes(ii).CorBad=Bad;
    
    
end
