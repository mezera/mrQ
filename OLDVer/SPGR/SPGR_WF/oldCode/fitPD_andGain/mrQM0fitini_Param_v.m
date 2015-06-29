function [x0,inDat,bm,szB,skip coils_list SDD]= mrQM0fitini_Param_v(opt,fb,M,coils,strongbiasLevel)
%
%function [x0,inDat,bm,szB,skip coils_list SDD]= mrQM0fitini_Param_v(opt,fb,M,coils,strongbiasLevel)
%
% INPUTS:
% opt - this is optmization structure that was passed from
% fb - ocation of the box in grid space. this is location of the box (this is not x,y,z location in image space but
% grid of boxes we made by mashgrid in mrQ_fitPD_multicoil
% M - the multi coils (4D) M0 image
% coils - howmany coils in the mesurament
% strongbiasLevel - a flag for steep coil gain like in the case of 32
% chanels
%
% OUTPUTS:
% x0 - intial guess for the parameters to fit
% inDat -the data to fit
% bm - brain mask for this region (box)
% szB - size opf the box
%skip - if this data need to be skip or not (1 yes)
%coils_list - the coil list we took from the M0
%SDD - variabilty within the data estimation (not in use)

%% I. Initialization
%get freesurfere segmentation
fileFS=fullfile(opt.outDir,['T1w_tissue.nii.gz']);
if(exist(fileFS,'file'))
    mask = readFileNifti(fileFS);
    mask=double(mask.data);
else
    disp(['error , can not find thoutcoilse file: ' fileFS]);
    error
end;
% get the kind of the problem
if(~exist('strongbiasLevel','var') || isempty(strongbiasLevel))
    if coils==32,
        strongbiasLevel=1;
        %32chancel have alot of bias and a lot of redundene information we try to deal with it. in the future more coils withthe same pattren should be added
    else
        strongbiasLevel=0;
    end;
end;


coils_list=0;

%get the location of the box we work on in image space (x,y,z)
[Xx Yy Zz,skip]=MrQPD_boxloc(opt,fb);

if skip==1;
    % somthing wrong no data no brain ...
    x0=0;known=0;inDat=0;bm=0;szB=0;
    return
end;


%% get free surfare mask brain mask and data in the box

% save the freesurfere segmentation in the box space
mask1=mask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));

% save brain mask in the box space
bm=opt.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
szB=size(bm);

%not in use but 32chancel might havve bliund spot that make every thing harder
%Bedge=length(find(bm1))./length(find(bm));
%Tres=squeeze(max(max(max((M.data(:,:,:,:))))))./60; %the is blind spots and low signal that we like to mask out in each coil

%get the M0 data in the box space
% use the data or a sqrt version of it (the sqrt is easyer to fit but we might loss
% details).
if opt.sqrtF==1
    box(:,:,:,:)=sqrt(double(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),:)));
else
    box(:,:,:,:)=(double(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),:)));
end
for i=1:coils;
    if opt.sqrtF==1
        
        tmp=box(:,:,:,i).^2;
    else
        tmp=box(:,:,:,i);
    end
    % if the signal is very low better not to use it (for sharp coils like
    % 32 head coils)
    tmp1=M.data(:,:,:,i);
    Tres(i)=prctile(tmp1(opt.brainMask),20);
    rate(i)=mean(tmp(bm));
    if find(tmp(bm)<Tres(i) & tmp(bm)>0 & strongbiasLevel==1 )
        
        % if find(tmp(bm1)<Tres(i) & tmp(bm1)>0 & strongbiasLevel==1 & Bedge==1)
        if (length(find(tmp(bm)<Tres(i) & tmp(bm)>0 ))./length(find(bm))>0.02 )
            % if (length(find(tmp(bm1)<Tres(i) & tmp(bm1)>0 ))./length(find(bm1))>0.02 )
            
            rate(i)=0;
        end;
    end
    
end;
%rate the coils by the there signal
[v In]=sort(rate,'descend');

%find the coils that we like to use
outcoils=min(coils,opt.numIn.*2);
outcoils=min(outcoils,length(find(rate>0)));

%clean out layer in the brain mask region in the box. we will take those
%point form the box and not fit accurding to those voxels
for i=1:outcoils,
    tmp=box(:,:,:,In(i));
    no=find(isnan(tmp));
    bm(no)=0;
    no=find(isinf(tmp));
    bm(no)=0;
    Mm=median(tmp(bm));
    Sd=std(tmp(bm));
    
    no=find(tmp>(Mm+3*Sd));
    bm(no)=0;
    no=find(tmp<(Mm-3*Sd));
    bm(no)=0;
    no=find(tmp<=0);
    bm(no)=0;
    
    
end;


%check that we still have voxel to work with
if (length(find(bm))<200 | outcoils==0) ;
    x0=0;known=0;inDat=0;bm=0;szB=0;skip=1;
    return
end;
for i=1:outcoils,
    tmp=box(:,:,:,In(i));
    inDat(:,i)=tmp(bm);
end;

%% intialize the fitting and check for bad coils

% the coils that are fine (now it still all of them)
posible=ones(outcoils,1);

for i=1:outcoils,
%go over coil and try to fit polynomial to the data. this is first
%estimation of the coil gain
    
%fit poylnomyals
    [params1,gains,rs ] = fit3dpolynomialmodel(box(:,:,:,In(i)),bm,opt.degrees);
    
% the calculate gain
    G=opt.Poly*params1';
  %in 3D
    G=reshape(G,size(bm));
    % the cacluate pd
    tmp=box(:,:,:,In(i))./G;
    %estimation of variabilty of the estimation
   
    stdC(i)=std(tmp(bm));
   
    ff(:,:,:,i)=tmp;
    if find(G(bm)<=0),
        %this is bad initial values better fit with out it
                posible(i)=0;
    end
% the intial estimation
    x0(i,:)=params1;
end;

% let check if the free surfarer can help us fit a better bias
%each tissue in the box
ti(1)=length(find(mask1==1));
ti(2)=length(find(mask1==2));
ti(3)=length(find(mask1==3));

most=(find(ti==max(ti))); %the tissue we have most in the box
if (ti(most(1))./length(find(bm))>0.03 && (max(ti)>40)); %is tissue mask is exsisted we will intiate the parameters with it
    
% go over the coils
    for i=1:outcoils,
        %fit poylnomyals
        [params1,gains,rs] = fit3dpolynomialmodel(box(:,:,:,In(i)),mask1==most(1),opt.degrees);
        % the calculate gain
        G=opt.Poly*params1';
          %in 3D
        G=reshape(G,size(bm));
        
        ff(:,:,:,i)=box(:,:,:,In(i))./G;
        if find(G(bm)<=0), %this is bad initial values better try to fit again with all the voxel
            % maskt=bm & logical(box(:,:,:,In(i)));
            [params1,gains,rs] = fit3dpolynomialmodel(box(:,:,:,In(i)),bm,opt.degrees);
            G=opt.Poly*params1';
            G=reshape(G,size(bm));
            ff(:,:,:,i)=box(:,:,:,In(i))./G;
            if find(G(bm)<=0),
                %this is bad initial values better fit with out it
                posible(i)=0;
            end
        end
        % the intial estimation
        x0(i,:)=params1;
        
    end
    
    
end


%make the box brain mask as in the same organization as the output data
bm=bm(:);
bm=repmat(bm,[1 outcoils]);



%% more check of the coils M0 data
%%

%this is the trechhold for coralation between coils
%if the coil are similar(coralte) we don't need it. as it will make the
%ability to define what is brain and what is coil harder
if strongbiasLevel==1;
    % a strong coralation level happan for small and close coils like the 32 chanel.
    %in this case we also have a lot of coils and we can try to avoid the
    cut=0.9;
else
    %where there are less coils we can't be to piky (8chancel coil)
    cut=0.95;
end

%1. check for M0 coil coralation in the box
% if the input are too corralate we don't won't them this is not adding
% information but may add biases.

%the coraltion
coefdat=tril(corrcoef(inDat),-1);
%loop to find to higly coralated and exclude one of them
if find(coefdat>cut);
    [C1(1,:) C2(1,:)]= ind2sub(size(coefdat),find(coefdat>0.9));
    
    nogood=find(posible==0);
    if length(nogood)~=0
        for i=1:length(nogood)
            
            C2(find(C2==nogood(i)))=0;
            C1(find(C1==nogood(i)))=0;
        end
    end
    
    
    while find(C1.*C2)
        tt1=[C1(find(C1.*C2)) C2(find(C1.*C2))];
        tt=sort(unique(tt1),'descend'); tt=tt(tt>0) ;
        if length(tt)>0
            clear count
            for i=1:length(tt)
                count(i)=length(find(tt1==tt(i)));
            end
            [d I]=max(count);
        end
        
        posible(tt(I))=0;
        C1(find(C1==tt(I)))=0;
        C2(find(C2==tt(I)))=0;
    end
end

%
%2. let check if this data is hard to fit. that can be
%becouse this the coil is very noise or can't be fited with the
%poyinoyals. This can happan in the coil blind spots where the data vary
%stiply or when the box is close to the coil (rare)
% the check is on the fitting of the poyinoyals indevitualy on each coil.
% we did that above when er initilazed the parameters. if this part is hard
% then it will be harder when we fit thogther with the other coils.
%we will also acludes look for nuns

% we avaluate the std of the fit small std mean easy to fir fit
Mstd=median(stdC(logical(posible)));

%no we go over coils and check
for i=1:outcoils;
 % no nuns please
    if (find(isnan(x0(i,:))))
        posible(i)=0;
    end;
    %no steep data we check the ratio between upper and lowest values
    if (min(inDat(:,i))/max(inDat(:,i))<0.01)
        posible(i)=0;
    end
    % no big std that is amarker of a dificulty to fit the box. (becouse of
    % stipness or nolise)
    if stdC(i)>2*Mstd;
        posible(i)=0;
    end
   
end

%clean the the coil we exclude for the ourput parameters
% any
if find(posible)
    bm=bm(:,find(posible));
    x0=x0(find(posible),:);
    inDat=inDat(:,find(posible));
    coils_list=In(find(posible));
else
    %if we have no coil left we will skip this box.
    skip=1;
end

