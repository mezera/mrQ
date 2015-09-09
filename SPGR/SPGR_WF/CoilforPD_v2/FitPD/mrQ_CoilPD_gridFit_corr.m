function mrQ_CoilPD_gridFit_corr(IDnum,jumpindex,jobindex)
% mrQ_CoilPD_gridFit(IDnum,jumpindex,jobindex)
%
%  This function is called by the SunGrid. It loads the relevant data and
%  fits the PD and coils bias of the M0 image region (volume, also called a
%  "box"). The box is a location (hundreds to thousands of voxels).  The
%  idea is to use the information from each coil. We try to find the PD
%  that is similar to all the coils' images, and fit the coil bias
%  that is different for each coil.
%
%  The fit is done in a few steps:
%  1. Get the box data.
%  2. Fit a subset of the coils information
%  3. Save the resuts.
%
%   ~INPUTS~
%             opt:   This is the optimization structure which was passed
%                         from mrQ_fitPD_multicoil. It contains all the
%                         needed information.
%       jumpindex:   How many boxes this grid call we fit (bookkeeping)
%        jobindex:   The number of boxes it will start when calling the
%                         grid (bookkeeping)
%
%  ~OUTPUTS~
%  This function will save an output file with fitted parameters in a tmp
%  directory. This will be used later by mrQfitPD_multiCoils_M0 to make the
%  PD map
%
% SEE ALSO:
% mrQ_PD_multicoil_RgXv_GridCall
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization

mrQpath= mrQ_getPath(IDnum);
load(mrQpath);
load(mrQ.opt_logname);

%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%check that this box contains brain data
if ed>length(opt.wh), ed=length(opt.wh);end;

nIteration=ed-st+1;
%initialize parameters and saved outputs

% Get the M0 and T1 information

%multi-coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;

% %T1
% T1=readFileNifti(opt.T1file);
% T1=T1.data;

%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

%seg mask
seg=readFileNifti(opt.segfile);
seg=seg.data;

smoothkernel=opt.smoothkernel;

% the poly basis to fit the coil gains
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);
maxCoil=opt.maxCoil;
minCoil=opt.minCoil;
useCoil=opt.useCoil;
nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

% initiate the saved parameters
fb=zeros(nIteration,1,3);
gEst=zeros(nPolyCoef,maxCoil,nIteration);
resnorm=zeros(nIteration,1);
exitflag=zeros(nIteration,1);
skip=zeros(nIteration,1);

Clists=zeros(maxCoil, nIteration);
Clists2=zeros(maxCoil, nIteration);
ResidErr=zeros(nIteration,1);

Iter=0;

options = optimset('Display','off',...  %'iter'final
    'MaxFunEvals',Inf,...
    'MaxIter',200,...
    'TolFun', 1e-6,...
    'TolX', 1e-6,...
    'Algorithm','levenberg-marquardt');

%%  II. Go over it, box by box

for ii= st:ed,
    %run over the box you like to fit
    clear M01  t1  BM1  SZ M0_v R1basis PDinit Segmask g0 G0 mask1 X_valdationErrF   X_gEstF  best2 Segmask
    Iter= Iter+1;
    tic
    %find the x,y,z location of the box (this is not x,y,z location in image space but
    %grid of boxes we made by meshgrid in  mrQ_PD_multicoil_RgXv_GridCall.m
    [fb(Iter,1,1), fb(Iter,1,2), fb(Iter,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, t1, BM1, SZ, skip(Iter), Segmask]= mrQ_GetM0_boxData(opt,[],M0,BM,fb(Iter,1,:),smoothkernel,seg);
    
    if  skip(Iter)==1
        %         disp(['skipping box ' num2str(ii) ' bad data'])
        
    else
        
        M0_v = reshape(M01, prod(SZ(1:3)), SZ(4));
%         coefdat = tril(corrcoef(M0_v(BM1,Clist)),-1); %% Clist not yet defined!
        
        %% III. Select the coils to fit
        %
        if length(useCoil)>SZ(4);
            useCoil=useCoil(1:SZ(4));
        end
        
        Clist=mrQ_select_coils(maxCoil,max(useCoil),M0_v(BM1,:));
        nCoils=length(Clist);
        Clists(1:nCoils,Iter)=Clist;
        
        % find an alternative coil list
        Mtmp=M0_v(BM1,:);
        Mtmp(:,Clist)=1;
        Clist2=mrQ_select_coils(maxCoil,max(useCoil),Mtmp);
        Clists2(1:nCoils,Iter)=Clist2;
        
        clear Mtmp;
        
        %% IV. Initiate the search parameters
        
        
        [ PDinit,g0] =Get_PDFit_InIt(3,M0_v(:,Clist),pBasis,[],Segmask)      ;
        
        %%
        %% moved this from lines 120
        coefdat = tril(corrcoef(M0_v(BM1,Clist)),-1);
        %%
        RegWeight  = 1000;
        XvalidationMask = logical(M0_v(BM1,Clist));
        
        [gEst(:,:,Iter), resnorm(Iter),~,exitflag(Iter)] = ...
            lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,double(M0_v(:,Clist)),pBasis,nCoils,RegWeight,BM1,double(coefdat),XvalidationMask),g0,[],[],options);
        
        [PDfit, Gn] = pdEstimate(M0_v(BM1,Clist),pBasis, gEst(:,:,Iter));
        %%        % get the  PD with a different set of coils
        
        coefdat = tril(corrcoef(M0_v(BM1,Clist2)),-1);
        
        [g, resnorm1,dd1,exitflag1] = ...
            lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,double(M0_v(:,Clist2)),pBasis,nCoils,RegWeight,BM1,double(coefdat),XvalidationMask),gEst(:,:,Iter),[],[],options);
        
        [PDfit2, Gn] = pdEstimate(M0_v(BM1,Clist2),pBasis, g);
        
        
        %%
        tmp=zeros(SZ(1:3));
        tmp(BM1)=PDfit2(:);
        PDfit2=tmp;
        tmp(BM1)=PDfit(:);
        PDfit=tmp;
        
        PDfit2  = reshape(PDfit2,SZ(1:3));
        PDfit  = reshape(PDfit,SZ(1:3));
        
        ErrMap=PDfit-PDfit2;
        [~,~,ResidErr(Iter)] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);
        
        toc
        %
    end
end

%% VI. Save
name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'gEst','resnorm','exitflag','st','ed','skip','fb' ,'Clists','Clists2','ResidErr')

