function mrQ_PDfit_ParallerCoils_Gridcall(optName,jobindex);
% mrQ_PDfit_ParallerCoils_Gridcall(opt,jobindex);
%INPUT
%optName %       a file name that has a stracture (opt) with information that is saved by
%                mrQ_fitPD_multicoil.m for the grid call
%jobindex       the Opt is the input for all the different grid calles. jobindex
%               indetify the job need to be run in this grid call

%outPut
% the function will save output file with the results
%
% look also at mrQ_PDfit_ParallerCoils_Call.m

load (optName);
% lets find what coil couple we are working on
A=size(opt{1}.couples);
B=(find(opt{1}.couples));
[i j]=ind2sub(A,B(jobindex)); 
j=opt{1}.couples(i,j);
% lets load the data
 M0=readFileNifti(opt{1}.dat);
in1=M0.data(:,:,:,i);
in2=M0.data(:,:,:,j);

in2cut1=M0.data(:,:,:,opt{1}.coil2coil(i));
in2cut2=M0.data(:,:,:,opt{1}.coil2coil(j));

clear M0;
BM=readFileNifti(opt{1}.BMfile);
BM=logical(BM.data);

%% lets find the mask of coil1
med=median(in1(BM));
up=prctile(in1(BM),opt{1}.prctileClip); %we clip the highest SNR voxels bexouse they are imposible to fit with polynomials
mask=BM;
mask(in1<med)=0; % we clip the noise part below of the median value (maybe this need to be different for differnt coils
mask(in1>up)=0;

% we will intersect it with a reference coil
med1=median(in2cut1(BM));
ratio=(in1./med)./(in2cut1./med1);
M0mask(:,:,:,1)=mask & ratio>0.3 ;

clear med up mask med1 ratio in2cut1


%% lets find the mask of coil1
med=median(in2(BM));
up=prctile(in2(BM),opt{1}.prctileClip); %we clip the highest SNR voxels bexouse they are imposible to fit with polynomials
mask=BM;
mask(in2<med)=0; % we clip the noise part below of the median value (maybe this need to be different for differnt coils
mask(in2>up)=0;

% we will intersect it with a reference coil
med1=median(in2cut2(BM));
ratio=(in2./med)./(in2cut2./med1);
M0mask(:,:,:,2)=mask & ratio>0.3 ;

clear med up mask med1 ratio in2cut2

[Poly,str] = constructpolynomialmatrix3d(size(BM),find(ones(size(BM))),opt{1}.degrees);
TM=readFileNifti(opt{1}.TM);
TM=TM.data;
%% map1
    mask=TM ==2 & in1;
    [params1,gains,rs] = fit3dpolynomialmodel(in1,mask,opt{1}.degrees);
    
    % we check if there are point that are off with this model (we
    % won't use it)
    G=Poly*params1';
    G=reshape(G,size(BM));
    % up date the mask and  exclude crazy points 50% more or less then WH %
    % we also exclude point that are end up having greater value after
    % corrections
    M0mask(:,:,:,1)=M0mask(:,:,:,1) & in1./G>0.5 & in1./G<1.5   & in1./G<in1;
    % re fit inital model
    mask=TM ==2 & M0mask(:,:,:,1);
    [params1,gains,rs] = fit3dpolynomialmodel(in1,mask,opt{1}.degrees);
    % up date the mask and  exclude crazy points ( we can do it more and
    % more  as a loop but maybe this is enghf.
        G1=Poly*params1';
            G1=reshape(G1,size(BM));

    M0mask(:,:,:,1)=M0mask(:,:,:,1) & in1./G1>0.5 & in1./G1<1.5  & in1./G1<in1;
    
    %save the inital parameters
    X0(1,:)=params1;
%% map1
    mask=TM ==2 & in2;
    [params1,gains,rs] = fit3dpolynomialmodel(in2,mask,opt{1}.degrees);
    
    % we check if there are point that are off with this model (we
    % won't use it)
    G=Poly*params1';
    G=reshape(G,size(BM));
    % up date the mask and  exclude crazy points 50% more or less then WH %
    % we also exclude point that are end up having greater value after
    % corrections
    M0mask(:,:,:,2)=M0mask(:,:,:,2) & in2./G>0.5 & in2./G<1.5   & in2./G<in2;
    % re fit inital model
    mask=TM ==2 & M0mask(:,:,:,2);
    [params1,gains,rs] = fit3dpolynomialmodel(in2,mask,opt{1}.degrees);
    % up date the mask and  exclude crazy points ( we can do it more and
    % more  as a loop but maybe this is enghf.
        G2=Poly*params1';
            G2=reshape(G2,size(BM));

    M0mask(:,:,:,2)=M0mask(:,:,:,2) & in2./G2>0.5 & in2./G2<1.5  & in2./G2<in2;
    
    %save the inital parameters
    X0(2,:)=params1;    
  
    
    %% join mask and data
    
mask= M0mask(:,:,:,1) & M0mask(:,:,:,2) & abs(1 - (in2./G2) ./ (in1./G1))<0.3;

M(:,1)=in1(mask);
M(:,2)=in2(mask);

mask=mask(:);
 % raw correlations
clear   G G1 G2 in1 in2 str BM M0mask 
MM=M;
MM(isnan(MM))=0;

    coefdat=corrcoef(MM);
    coefdat=coefdat(1,2);
    
clear MM

%% Lets Fit

a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'iter' ,'Tolx',1e-6,'MaxIter',200,'MaxFunEvals',inf);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'on','Tolx',1e-6,'MaxIter',200,'MaxFunEvals',inf);%'TolF',1e-12
    
end
name=[ opt{1}.name '_' num2str(jobindex) '_']


[res, resnorm,dd1,exitflag] = lsqnonlin(@(par) mrQ_errParallGainFitCouple(par,double(M),Poly,double(coefdat),mask),double(X0),[],[],options);
 

%name=[ opt{1}.name '_' jobindex '_'];

%check where the two fit are good both values are similar
%make a waithed sum PD map;
%output(:,1)=location;
%output(:,2)=values (PDmap);
% save( degree,output,res1)
% i still need to think about what i'm saving:

%% save the output 1 mask 2 first coil PD 3 secound coil PD
output(1,:)=find(mask);

G=Poly*res(1,:)';
output(2,:)=M(:,1)./G(mask);

G=Poly*res(2,:)';
output(3,:)=M(:,2)./G(mask);

save(name,'output','resnorm','exitflag','i','j')


