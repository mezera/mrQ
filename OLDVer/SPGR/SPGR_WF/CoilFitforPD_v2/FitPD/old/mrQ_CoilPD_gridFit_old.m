function mrQ_CoilPD_gridFit_old(opt,jumpindex,jobindex)
%
% mrQ_CoilPD_gridFit(opt,jumpindex,jobindex)
%  this function call by the sun grid it load the relavant data and fit the
%  PD and coils bias of M0 mage rigion (voulume).
% the imaging voulume region also call here "box". The box is a location (few voxel 100's to
% 1000's).  the idea is to use the information in each coils. we tery to find
% to PD that is similar to all the coils images and fit the coil bias that
% is diferent for each coil.
%  The fit is done in few steps.
%  1. we get the box data.
% 2. we fit subset of coils informatio nwith regularizartion (of T1) and crossvalidation to it's quallety
% 4. the best cross validation regularization is used to fit all the data.
%4. we save the resuts
%
% INPUTS:
%       opt - this is optmization structure that was passed from
%       mrQ_fitPD_multicoil and it have all the needed information
%       jumpindex - how many boxes this grid call we fit (book keeping)
%       jobindex  - the number of box it will start whe ncalling the grid
%       (book keeping)
%
% OUTPUTS:
%  save an output file with fitted parameters in a tmp directorry
%   this will be used lster by mrQfitPD_multiCoils_M0 to make the PD map

% SEE ALSO:
% mrQ_PD_multicoil_RgXv_GridCall
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization



%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%check that this box have brain data
if ed>length(opt.wh), ed=length(opt.wh);end;

nIteration=ed-st+1;
%intilazie parameters and saved outputs

% Get the M0 and T1 information

%multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
%T1
T1=readFileNifti(opt.T1file);
T1=T1.data;


%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

%seg mask
seg=readFileNifti(opt.segfile);
seg=seg.data;


smoothkernel=opt.smoothkernel;


% thepoly basis to fit the coil gains
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);
maxCoil=opt.maxCoil;
minCoil=opt.minCoil;
useCoil=opt.useCoil;
nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

% initite the saved parameters
fb=zeros(nIteration,1,3);
gEst=zeros(nPolyCoef,maxCoil,nIteration);
resnorm=zeros(nIteration,1);
exitflag=zeros(nIteration,1);
skip=zeros(nIteration,1);
X_valdationErrSN=zeros(nIteration,1);
X_valdationErr=zeros(2,length(opt.lambda),nIteration);
BestReg=zeros(nIteration,2);
Clists=zeros(maxCoil, nIteration);
Iter=0;

%%  II. go over the box the boxs


for ii= st:ed,
    %run over the box you like to fit
    clear M01  t1  BM1  SZ M0_v R1basis PDinit Segmask g0 G0 mask1 X_valdationErrF   X_gEstF  best2 Segmask
    Iter= Iter+1;
    tic
    %find the x,y,z location of the box (this is not x,y,z location in image space but
    %grid of boxes we made by mashgrid in  mrQ_PD_multicoil_RgXv_GridCall.m
    [fb(Iter,1,1), fb(Iter,1,2), fb(Iter,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, t1, BM1, SZ, skip(Iter), Segmask]= mrQ_GetM0_boxData(opt,T1,M0,BM,fb(Iter,1,:),smoothkernel,seg);
    
    if  skip(Iter)==1
        disp(['skipping box ' num2str(ii) ' bad data'])
        
    else
        
        M0_v = reshape(M01, prod(SZ(1:3)), SZ(4));
        %make a R1 regularazation matrix
        R1basis(1:nVoxels,1) = 1; R1=1./(t1(:)); R1basis(:,2) =R1;
        R1basis=double(R1basis);
        PDinit=1./(R1*0.42+0.95); %this is the teortical T1 PD relationship see reference at Mezer et. al 2013
        
        
        %% Select the coils to fit
        %
        if length(useCoil)>SZ(4); useCoil=useCoil(1:SZ(4));end
        c=nchoosek(useCoil,maxCoil);
        Cor=ones(maxCoil,size(c,1))*100;
        for k=minCoil:maxCoil
            
            c=nchoosek(useCoil,k);
            
            for kk=1:size(c,1)
                A=(corrcoef(M0_v(BM1,c(kk,:))));
                %    Cor(k,kk)=(sum(A(:))-k)/2;
                Cor(k,kk)=(sum(abs(A(:)))-k)/((length(A(:))-k)*0.5);
                
                Cloc(k,kk,1:k)=c(kk,:);
            end
        end
        
        
        [v ind]=sort(Cor(:)); %sort the coils by minimum corralation
        [xx yy]=ind2sub(size(Cor),ind(1)); % find the combination with minimal corralation
        Clist=[squeeze(Cloc(xx,yy,1:xx))']; % get the coil list
        nCoils=length(Clist);
        Clists(1:nCoils,Iter)=Clist;
        %% Let's make sure we have no strange voxel
        
        for jj=1:nCoils
            M = M0_v(:,Clist(jj)) ; % Raw estimate
        Bad=isnan(M) | isinf(M) | M==0;
        BM1(Bad)=0;
        end
         Bad=isnan(R1) | isinf(R1) | R1==0;
        BM1(Bad)=0;
        

        %% intiate the search parameters
        
        
        
        G  = zeros(nVoxels,nCoils);
        g0 = zeros(nPolyCoef,nCoils);
        
        mask1 =find(BM1);
        
        for jj=1:nCoils
            G(mask1,jj)  = M0_v(mask1,Clist(jj)) ./ PDinit(mask1);         % Raw estimate
            g0(:,jj) = pBasis(mask1,:) \ G(mask1,jj);  % Polynomial approximation
        end
        
        if any(isnan(g0(:)));        g0(isnan(g0))=0;        end
        %%  X-Validation Fit of coil Gain with T1 regularization
        
        [X_valdationErrF,  X_gEstF]=pdX_valdationLoop_2(opt.lambda,3,M0_v(BM1,Clist), pBasis(BM1,:),R1basis(BM1,:),g0,Segmask(BM1));
        
        %       [X_valdationErrF,  X_gEstF]=pdX_valdationLoop_2(opt.lambda,opt.Kfold,M0_v(BM1,opt.Clist), pBasis(BM1,:),R1basis(BM1,:),g0,Segmask(BM1));
        X_valdationErr(:,:,Iter)=X_valdationErrF;
        % minimum X_valdation error
        best1 = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:))) % sum of abs err
        best2 = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)))% RMSE
        %figure;  semilogy(opt.lambda,X_valdationErrF(2,:),'*-'); X_valdationErrF(2,:)./min(X_valdationErrF(2,:))
        %figure;  semilogy(opt.lambda,X_valdationErrF(1,:),'*-'); X_valdationErrF(1,:)./min(X_valdationErrF(1,:))
        BestReg(Iter,:)=[best1 best2];
        % i cases we get no meaningfull Xvalidation with the
        % regularization we should be worry!!!
        % if there is not enght data we many X_valdation Err are
        % almost similar we can cause to take the begest one. this
        % need to be implamented
        if min(X_valdationErrF(2,:))*2 > X_valdationErrF(2,end)
            X_valdationErrSN(Iter)=1;
            disp(['no clear  X-validation improvment with regularization'])
        end
        
        %  [PDfit,RMSE1]=pdCoilSearch_T1reg( lambda1(best1),M0_v(:,Clist),pBasis,R1basis,X_gEstF(:,:,1,best),[],[],PDsim);
        [PDfit,~,G,gEst(:,1:xx,Iter), resnorm(Iter),exitflag(Iter) ]=pdCoilSearch_T1reg( opt.lambda(best2),M0_v(BM1,Clist),pBasis(BM1,:),R1basis(BM1,:),X_gEstF(:,:,1,best2),Segmask(BM1));
        
        
        toc
        %
    end
end

name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'gEst','resnorm','exitflag','st','ed','skip','fb' ,'X_valdationErr','X_valdationErrSN','BestReg','Clists')


