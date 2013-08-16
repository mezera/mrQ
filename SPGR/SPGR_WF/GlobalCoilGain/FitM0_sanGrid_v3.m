function FitM0_sanGrid_v3(optName,jobindex)
%
% FitM0_sanGrid_v2(opt,jumpindex,jobindex)
%  this function call by the sun grid it load the relavant data and fit the
%  PD and coils bias of image rigion.
%the image region also call here box. is a location (few voxel 100's to
%1000's).  the idea is to use the information in man coils. we tery to find
%to PD that is similar to all the coils images and fit the coil bias that
%is diferent for each coil.
%  The fit is done in few step. 1. we get the box data.
%2. we fit subset of coils (nested function )
%3. we do qualtiy check to see if it worked o.k.  
%4. we save the resutl
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
% mrQ_fitPD_multicoil
% (C) Stanford University, VISTA
%
%

%% I. Initialization
load (optName);
jumpindex=opt{1}.jumpindex;
% use only one cpu or it will be very hard to work on the computer we used.
%maxNumCompThreads(1);
%maxnumcompthreads(1);

%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%cheack that this box have brain data
if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;

%intilazie parameters and saved outputs
nn=1;
res=[];resnorm=0;
exitflag=-100;
coilList=0;skip=0;

%%
%load the multi coil voulume
M0=readFileNifti(opt{1}.dat);
coils=size(M0.data,4);

coilList=zeros(coils,ed-st+1);
badbox=zeros(ed-st+1,1);

prctileClip=opt{1}.prctileClip;
%% load  the data

BMfile=opt{1}.HMfile;

BM=readFileNifti(BMfile);BM=logical(BM.data);
%book keeping:
%number of coils mesured

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
clear in1 in2 med med1 ratio
    nn=0;
T1=readFileNifti(opt{1}.T1file); T1=T1.data;
degrees=opt{1}.degrees;
out=[];
for i= st:ed,
    nn=nn+1;
     %run over the box you like to fit
  
    tic
    %find the x,y,z location of the box (this is not x,y,z location in image space but
    %grid of boxes we made by mashgrid in  mrQ_fitPD_multicoil
    [fb(i,1,1) fb(i,1,2) fb(i,1,3)]=ind2sub(size(opt{1}.X),opt{1}.wh(i));
    
    %get the location of the box we work on in image space (x,y,z)
[Xx Yy Zz,skip]=MrQPD_boxloc(opt{1},fb(i,:));

if skip==~1
    T1box=T1(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
t1mask=T1box>0.4 & T1box<5;
    mask=ones(Xx(2)-Xx(1)+1,Yy(2)-Yy(1)+1,Zz(2)-Zz(1)+1);
    
    tvx=cumprod(size(mask)); tvx=tvx(end);
    [Poly,str] = constructpolynomialmatrix3d(size(mask),find(ones(size(mask))),degrees);
   
    clear X0 in Boxmask
    j=0;
    for ii=1:coils
        inmask=M0mask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),ii) & t1mask;
        
        if length(find(inmask))/tvx>0.25
                        j=j+1;

            in=M0.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),ii);
            
            [params1,gains,rs] = fit3dpolynomialmodel(in,inmask,degrees);
            G=Poly*params1';
            G=reshape(G,size(mask));
            % up date the mask and  exclude crazy points 50% more or less then WH %
            % we also exclude point that are end up having greater value after
            % corrections
            Boxmask(:,:,:,j)=inmask & in./G>0.5 & in./G<1.5   & in./G<in;
            % re fit inital model
            inmask= Boxmask(:,:,:,j);
            [params1,gains,rs] = fit3dpolynomialmodel(in,inmask,degrees);
            
            % up date the mask and  exclude crazy points ( we can do it more and
            % more  as a loop but maybe this is enghf.
            G=Poly*params1';
            G=reshape(G,size(mask));
            Boxmask(:,:,:,j)=inmask & in./G>0.5 & in./G<1.5   & in./G<in;
            
            X0(:,j)=params1';
            coilList(ii,nn)=1;
        end
    end
    % let check the  coils ovelap mask 
    if j>5
    mask=(sum(Boxmask,4));
    mask(mask<5)=0;
    mask=logical(mask>5);
    else
        mask=1;
    end
    if ~length(find(mask))/tvx>0.25
        skip=1;
    end
end; 



%% fit or skip the box
if skip==1 ||  j<6  % skip
        badbox(nn)=1
        x0=0;known=0;inDat=0;bm=0;szB=0;
        
else %fit
      % we will use stocastic lsq fit
          [Poly,str] = constructpolynomialmatrix3d(size(mask),find((mask)),degrees);

  
 M=nan(length(find(mask)),j);
CL=find(coilList(:,nn));
for ii=1:j;
    
    inMask= Boxmask(:,:,:,ii);
    in=M0.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),CL(ii));
            
    in(~inMask)=nan;
    M(:,ii)=in(mask);
end

med=median(T1box(mask));
T1box=T1box(mask);
mask2=T1box>med*.9 & T1box<med*1.1;

raw=M;
        bad = isnan(raw);
        raw(bad)=0;
        coefdat =  tril(corrcoef(raw),-1);
      
        raw(~mask2,:)=0;

              coefdatT1 =  tril(corrcoef(raw),-1);

clear raw bad  mask

%
     a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'iter' ,'Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'iter','Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);%'TolF',1e-12
    
end



[res, resnorm,dd1,exitflag] = lsqnonlin(@(par) errParallGainFit_v2(par,M,Poly,coefdat,coefdatT1,mask2),double(X0),[],[],options);

       
out(nn).res=res;
out(nn).resnorm=resnorm;
out(nn).exitflag=exitflag;
out(nn).CL=CL;

    end
    
end
name=[ opt{1}.name '_' num2str(st) '_' num2str(ed)];
save(name,'out','badbox','st','ed', 'str')










