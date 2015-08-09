function mrQ_CoilPD_gridFit_PD_T1(opt,jumpindex,jobindex)
%
% mrQ_CoilPD_gridFit_PD_T1(opt,jumpindex,jobindex)

%
% INPUTS:
%         opt:   This is the optimization structure that was passed from
%                   mrQ_fitPD_multicoil. It has all the needed information.
%   jumpindex:   How many boxes we fit in this grid call (bookkeeping)
%    jobindex:   The number of boxes it will start when calling the grid
%                  (bookkeeping)
%
% OUTPUTS:
%                The function will save an output file with fitted
%                parameters in a tmp directory. This will be used later by
%                mrQfitPD_multiCoils_M0 to make the PD map.

% SEE ALSO:
% mrQ_PD_multicoil_RgXv_GridCall
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization

%Find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%Check that this box has brain data
if ed>length(opt.wh), ed=length(opt.wh);end;

nIteration=ed-st+1;
%Initialize the parameters and saved outputs

% Get the M0 and T1 information

% Multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;

%T1
T1=readFileNifti(opt.T1file);
T1=T1.data;

%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

BM(T1>3000)=0; % Clear areas that are not GM or WM.

%seg mask
seg=readFileNifti(opt.segfile);
seg=seg.data;

smoothkernel=opt.smoothkernel;

% The poly basis to fit the coil gains
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);

nVoxels=size(pBasis,1);
nPolyCoef=size(pBasis,2);

% Initiate the saved parameters
fb=zeros(nIteration,1,3);
gEst=zeros(nPolyCoef,nIteration);
resnorm=zeros(nIteration,1);
exitflag=zeros(nIteration,1);
skip=zeros(nIteration,1);

Iter=0;

%%  II. Go over it, box by box

for jj= st:ed,
    %run over the box you like to fit
    clear M01  t1  BM1  SZ M0_v R1basis PDinit Segmask g0 G0 mask1
    Iter= Iter+1;
    tic
    %Find the x,y,z location of the box.
    % (This is not the x,y,z location in image space, but rather the grid
    % of boxes we made by meshgrid in  mrQ_PD_multicoil_RgXv_GridCall.m)
    [fb(Iter,1,1), fb(Iter,1,2), fb(Iter,1,3)]=ind2sub(size(opt.X),opt.wh(jj));
    
    % Get all the relevant box data for the fit
    [M01, t1, BM1, SZ, skip(Iter), Segmask]= mrQ_GetM0_boxData(opt,T1,M0,BM,fb(Iter,1,:),smoothkernel,seg);
    M0_v = M01(:);
    t1=t1.*1000;
    R1=1./(t1(:));
    
    
    if  skip(Iter)==1
       % disp(['skipping box ' num2str(jj) ' bad data'])
        
    else
 
        %% Fit
        %% Initiate the multi-box fit parameters  1/PD=A +B/T1
        A=zeros(1,7);B=zeros(1,7);
        A(1) = 0.916 ; B(1) = 436; %literature values
        
        % a basis for estimating the new A and B.
        R1basis(:,2)=1./t1(BM1);
        R1basis(:,1)=1;
        
        %iterate to converge A and B
        for ii=2:7
            PDp=1./(A(ii-1)+B(ii-1)./t1);
                      %  PDp=PDp./median(PDp(BM1)); scale

            % the sensitivity the receive profile
            RPp=M01./PDp;
            
            % Raw estimate
            g = pBasis(BM1,:) \ RPp(BM1);  % Polynomial approximation
            RPi=pBasis*g;
            
            % calculate PD from M0 and RP
            PDi=M01(:)./RPi;
            
            % solve for A B given the new PD estimation
            % ( 1./PDi(BM1) )= A* R1basis(:,1) + B*R1basis(:,2);
            
            co     = R1basis \ ( 1./PDi(BM1) );
            A(ii)=co(1);
            B(ii)=co(2);
            
        end
        gEst(:,Iter)=g;
    end
end

name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'gEst','st','ed','skip','fb')

