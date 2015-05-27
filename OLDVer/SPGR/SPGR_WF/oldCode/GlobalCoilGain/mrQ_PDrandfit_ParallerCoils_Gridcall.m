function mrQ_PDrandfit_ParallerCoils_Gridcall(optName,jobindex)
% mrQ_PDfit_ParallerCoils_Gridcall(opt,jobindex);
%INPUT
%optName %       a file name that has a stracture (opt) with information that is saved by
%                mrQ_fitPD_multicoil.m for the grid call
%jobindex       the Opt is the input for all the different grid calles. jobindex
%               indetify the job need to be run in this grid call

%outPut
%



load (optName);

degrees=opt{1}.degrees ;
outDir=opt{1}.outDir;


BMfile=opt{1}.BMfile;

M0cfile=opt{1}.M0cfile;

prctileClip=opt{1}.prctileClip;
%% load  the data

M0=readFileNifti(M0cfile);
BM=readFileNifti(BMfile);BM=logical(BM.data);
%book keeping:
%number of coils mesured
coils=size(M0.data,4);
%number of coil we will use
coilN=min(opt{1}.Ncoils,coils);

if notDefined('prctileClip')
    if coils>10; % i just guess here this need to be check with the spesipic coil in use
        prctileClip=98;
    else
        prctileClip=99;
    end
end


%% mask the data from crazy values
for i=1:coils;
    
    in=M0.data(:,:,:,i);
    up=prctile(in(BM),prctileClip); %we clip the highest SNR voxels bexouse they are imposible to fit with polynomials
    med=median(in(BM));
    mask=BM;
    mask(in<med)=0; % we clip the noise part below of the median value (maybe this need to be different for differnt coils
    mask(in>up)=0;
    M0mask(:,:,:,i)=mask;
end
clear in up med mask

for i=1:coils;
    
    for j=1:coils
        in= M0mask(:,:,:,i)+M0mask(:,:,:,j);
        Val(j)=length(find(in==2));
    end
    Val(i)=0;
    coil2coil(i)=find(max(Val)==Val);
    in=M0.data(:,:,:,i);
    med=median(in(BM));
    in1=M0.data(:,:,:,coil2coil(i));
    med1=median(in1(BM));
    ratio=(in./med)./(in1./med1);
    M0mask(:,:,:,i)=M0mask(:,:,:,i) & ratio>0.3 ; %we clip the regions that a coil is very differnt (smaller) from the most smiliar coil to it. becouse this is a marker for area we can't fit (most of the time this is blind spot of the coil)
end



TM=readFileNifti(opt{1}.TM);
TM=TM.data;

%% one more masking by fit to White matter
[Poly,str] = constructpolynomialmatrix3d(size(BM),find(ones(size(BM))),degrees);

for i=1:coils;
    mask=TM ==2 & M0mask(:,:,:,i);
    [params1,gains,rs] = fit3dpolynomialmodel(M0.data(:,:,:,i),mask,degrees);
    
    % we check if there are point that are off with this model (we
    % won't use it)
    G=Poly*params1';
    G=reshape(G,size(BM));
    % up date the mask and  exclude crazy points 50% more or less then WH %
    % we also exclude point that are end up having greater value after
    % corrections
    M0mask(:,:,:,i)=M0mask(:,:,:,i) & M0.data(:,:,:,i)./G>0.5 & M0.data(:,:,:,i)./G<1.5   & M0.data(:,:,:,i)./G<M0.data(:,:,:,i);
    % re fit inital model
    mask=TM ==2 & M0mask(:,:,:,i);
    [params1,gains,rs] = fit3dpolynomialmodel(M0.data(:,:,:,i),mask,degrees);
    % up date the mask and  exclude crazy points ( we can do it more and
    % more  as a loop but maybe this is enghf.
    G=Poly*params1';
    G=reshape(G,size(BM));
    
    M0mask(:,:,:,i)=M0mask(:,:,:,i) & M0.data(:,:,:,i)./G>0.5 & M0.data(:,:,:,i)./G<1.5  & M0.data(:,:,:,i)./G<M0.data(:,:,:,i);
    
end;

clear G  params1 gains rs mask TM Poly str
%%  this is our all brain mask
mask=(sum(M0mask,4));
mask(mask<4)=0;
mask=logical(mask>4);
%% save only the relevant data
M=nan(length(find(mask)),coils);

for i=1:coils;
    inMask=M0mask(:,:,:,i);
    in=M0.data(:,:,:,i);
    in(~inMask)=nan;
    M(:,i)=in(mask);
end


clear M0 M0mask inMask in;
mask=mask(:);


%%


% make a randum roi of voxels in the image
for i=1:jobindex  % we like to seed it diferntily for each grid call (other wise the all ave the same "random set")
vox=randperm(length(M));
C=randperm(coils);
end
C=C(1:coilN); % the coils we will use (the one that have the most full representation of the random ROI we define

Nvalus=opt{1}.Nvalus;
M1=M(vox(1:Nvalus),:); % the subset of data that we will fit.

% find the voxel that have values in this roi
%[X,ID]=sort(sum(~isnan(M1),1),'descend');
X=sum(~isnan(M1),1);
toLow=find(X<Nvalus/3); %not enght voxel to fit
% we will check that we don't have coil with not mucha data
checkLow=1;
while checkLow
if intersect(C,toLow)
    C=randperm(coils);
C=C(1:coilN); % the coils we will use (the one that have the most full representation of the random ROI we define

else
  checkLow=0;
end
end

clear X toLow

%we making a list of the ratio of coils 1&2 1&3 ...9&10
Val=[];
for i=1:coilN-1
    clear V
    V(1,1:coilN-i)=i;
    V(2,1:coilN-i)=i+1:coilN;
    
    Val=[Val V];
end
In1=Val(1,:);
In2=Val(2,:);



% devide one coil by anther let us fit the coil bias ratio. when we devid the  PD
% component is divided out.

Ratio=M1(:,C(Val(1,:)))./M1(:,C(Val(2,:)));
clear Val
% we will mask for values that are relevant and not to big as well
mask2=abs(Ratio)<10 & ~isnan(Ratio) & ~isinf(Ratio);

WHM=find(mask);
% make the polynomyals we will fit to the ROI
[Poly,str] = constructpolynomialmatrix3d(size(BM),WHM(vox(1:Nvalus)),degrees);

% make a starting guess for the parameters by fit a polinomyar to the each coil in the ROI
for i=1:coilN
    in=zeros(size(BM));
    in(mask)=M(:,C(i));
    inMask=zeros(size(BM));
    inMask(WHM(vox(1:Nvalus)))=1;
    inMask=inMask & in>0;
    [params1,gains,rs] = fit3dpolynomialmodel(in,inMask,degrees);
    X0(:,i)=params1';
end
clear in imMask params1 gains rs WHM 


%% find the data correlation with in the ROI 

MM=M1(:,C);
MM(isnan(MM))=0;

coefdat=tril(corrcoef(MM),-1);

Mmask=MM>0;
clear MM M M1 mask

%% Lets Fit

a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'iter' ,'Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'iter','Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);%'TolF',1e-12
    
end



[res, resnorm,dd1,exitflag] = lsqnonlin(@(par) errcoilsRatio(par,Ratio,Poly,coefdat,mask2,In1,In2,Mmask),double(X0),[],[],options);


%% lets save
name=[ opt{1}.name '_' num2str(jobindex) '_'];

save(name,'res','str','exitflag','resnorm','C')

%% lets look at the results
%   [Poly1,str] = constructpolynomialmatrix3d(size(BM),find(ones(size(BM))),2);
% 
%     Gain=Poly1*res;
%     PD=nan(size(M));
% 
% 
%      for i=1:10
%           PD(:,i)=M(:,C(i))./Gain(mask,i);
%      end
%      PD1=nanmedian(PD,2);
% 
%          in=zeros(size(BM));
%          in(mask)=PD1;
%        PD1=nanstd(PD,[],2);
%      in1=zeros(size(BM));
%        in1(mask)=PD1;
%        clear PD1 
%

