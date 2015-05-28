function [B1,AnalysisInfo]=mrQ_smmothL_B1(B1epifile,AlignFile,outDir,B1file,xform,degree,AnalysisInfo,t1fileHM,SET1Fitfile,B1epiResidfile);
%
% [B1,AnalysisInfo]=mrQ_smmothL_B1(B1epifile,AlignFile,outDir,B1file,xform,degree,matlabpoolF,AnalysisInfo,t1fileHM,SET1Fitfile,B1epiResidfile);
%
% Perform smoth the SEIR B1 map and register to SPGR space
%
% INPUTS:
%    B1epifile         -     The B1 map path 
%    AlignFile         -     the path for the raw data that was used to the B1 Fit
%       outDir         -     Ouput directory where the resulting nifti files will
%                            be saved. 
%       xform          -     Transform
%       B1file         -     output B1 name
%       degree         -     polynomiyal degree for global fit of the smooth B1
%                           (defult 3)
%       AnalysisInfo   -      the fit  and registration information structure that will be used for
%                           registration and then it will be updated 
%       t1fileHM       -     a SPGR T1 map for registration  
%        SET1Fitfile   -    a SEIR T1 map for registration 
%       B1epiResidfile -    a residual file of the B1 fit that is used to find outlayers B1
%                       values
% OUTPUTS
% b1 -smooth B1 in SPGR space
%AnalysisInfo -info structure updated after the registration



%we smoth the SEIR B1 map by local and then feel the gap by global
%regrigions and then register back the smooth B1 to the SPGR space by Ants
%softwere

% See Also:
% mrQfit_T1M0_ver2


if (~exist('degree','var')|| isempty(degree)),
    degree=3;
end;


%%  B1 mask
% mask the B1 value that we don't like to use 
%1. more then 50% T1 variation between SEIR and SPGR is probabaly miss registration
%2. B1 in the edge of the serech space is point that didn't converge in the B1
%fit
%3.voxel with big fit residual for B1 or SEIR T1 are also out

%load the B1 fit
B1=readFileNifti(B1epifile);
SE_Xform=B1.qto_xyz;
pixdim=B1.pixdim ;
B1=double(B1.data);


%load the data that was us to the B1 fit
load (AlignFile)


%load the B1 fit residual
B1fitResid=readFileNifti(B1epiResidfile);
B1fitResid= B1fitResid.data;

%load the T1 fit residuals

load (SET1Fitfile);
SEIRResid=ll_T1(:,:,:,4); %fit residual
a_fitparam=ll_T1(:,:,:,2); %tissue parameter form theSEIR eqation
b_fitparam=ll_T1(:,:,:,3);%tissue parameter form theSEIR eqation
clear ll_T1

%no fot to big or small B1 (50%)
tisuuemask=  B1<1.5 & B1>.5 & ~isinf(SEIRResid) ;

% no for big residual
tt=~isinf(SEIRResid) & SEIRResid>0;

tisuuemask=tisuuemask & SEIRResid<prctile(SEIRResid(tt),98) & B1fitResid<prctile(B1fitResid(find(B1fitResid)),97) ...
    & Res{1}.im<5000 & Res{1}.im>350 & b_fitparam>prctile(b_fitparam(tt),1) &a_fitparam<prctile(a_fitparam(tt),99);

tisuuemask=tisuuemask & B1>prctile(B1(tisuuemask),1) & B1<prctile(B1(tisuuemask),99);
tisuuemask=tisuuemask & (Res{1}.im./(Res{2}.im*1000))>  0.3;

tisuuemask=tisuuemask & (Res{1}.im./(Res{2}.im*1000))<  3;




%% fit local regresiions

tmp1=zeros(size(tisuuemask));


sz=size(tisuuemask);
tt=ones(sz); 




%% II fit local regresiions
% we will smooth the B1 map by local regresiions
%the fit is done in tree steps
%1.we first estimate the voxel that have can be estimate with great confidance
%(<45% ) of the area under the filter
%local information 

%2. we fit  the others that are posible but with less confidance with
%more smooth (biger filter

%3. the one that are out of reach for local regration will be fitted by
%global polynomyial along the full B1 space (like a fillter along the all B1 map)

tmp1=zeros(size(tisuuemask));


sz=size(tisuuemask);
tt=ones(sz); 

%%%
%1. we check the local information by compareing to covariance of
%the filter with a all ones imgage (full information).

area=0.45;
%filter size
FS=30;
filter1=FS./pixdim;
[f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);

%define the avilable caverage
C1 = convn(tisuuemask,f1,'same');

%define the maxsimal caverage
CC1=convn(tt,f1,'same');

%the voxel that we will use
tisuuemask1=C1>max(CC1(:)).*area;

%where there are B1 estimation (x,y,z location)
[x y z]=ind2sub(size(tisuuemask),find(tisuuemask));

%were we will find the smooth B1  estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));

%local regrision
w1 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter1,[]);

%save the result
tmp1(find(tisuuemask1))=w1;




%%%
    %2. we fit all the others
%we will increase the filter size and reduce the area that is needed to be
%included


area=0.15;
%filter
FS=60;
filter1=FS./pixdim;
[f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);

%define the avilable caverage
C1 = convn(tisuuemask,f1,'same');
%where there are B1 estimation

%define the maxsimal caverage
CC1=convn(tt,f1,'same');

%the voxel that we will use
tisuuemask2=C1>max(CC1(:)).*area & tisuuemask1==0;

%where there are B1 estimation (x,y,z location)
[x y z]=ind2sub(size(tisuuemask),find(tisuuemask));

%were we will find the smooth B1 estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask2));

%local regrision
w2 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter1,[]);

%save the result
tmp1(find(tisuuemask2))=w2;




%%% 3.   we can do global polynomial if we can't fit localy some locations    
Imsz1=size(tmp1);

%make the 3D polynials
[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),degree);

%fit the polynials coefitents to the smooth B1 map
[params,gains,rs] = fit3dpolynomialmodel(tmp1,(tmp1>0),degree);

%reshape from a vector to 3D map
B1match = reshape(Poly1*params(:),Imsz1);

%find where the B1 value are missing
mask=logical(tmp1>0);

%fill the holls
tmp1(~mask)=B1match(~mask);


%%% last check
%we might get  scaling effect when the smoothing is not perefect , scale it
%back so the mean of the smooth map and the original is the same
Cal=median(B1(tisuuemask)./tmp1(tisuuemask));
tmp1=tmp1.*Cal;



%%% save the smooth B1 map
B1epifile1=fullfile(outDir,['B1_fit_lreg_N3.nii.gz']);
dtiWriteNiftiWrapper(single(tmp1), SE_Xform, B1epifile1);


%% III take B1  to the SPGR space and feel the gaps with global fits

%% we need to register the B1 map from the SEIR space to SPGR we are using ANTS softwere
 disp('ATNS non linear registration and Warp may take about 20 min ...')
 disp('...')

 %% run ANTS non linear parameter to register SPGR to epi
% [AnalysisInfo]=mrQ_NLANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1epifile1)
 [AnalysisInfo]=mrQ_RB_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1epifile1);

 %load the result B1
 B1_h=readFileNifti(AnalysisInfo.RB_B1_epi_spgr);
 B1_h=B1_h.data;

 %%%% last things before this B1 is ready
 % feel any and funny values that were added in the registration
% we will use again the global polynomial fit in the SPGR space
 
%%



%let feel the gaps in the linear regration done in epi space that was Warp
%to SPGR with polynomial model of the values. 



%the first slice is always empty so what ever the value in there it need to be change.
%ants some time feel empty imafe spaces with some value that is randome but
%constant

out=mean(mean(B1_h(:,:,1)));
if isnan(out);
    out=0;
end
%find when the B1 is defently worng
mask=B1_h>0.3 & B1_h<1.7 & ( B1_h<out.*0.9999 | B1_h>out.*1.0001) ;

%polynomial fit of the values. for the missing values

Imsz1=size(B1_h);

%make the 3D polynials
[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),3);

%fit the polynials coefitents to the smooth B1 map
[params,gains,rs] = fit3dpolynomialmodel(B1_h,mask,3);

%reshape from a vector to 3D map
B1match = reshape(Poly1*params(:),Imsz1);

% feeling the gaps and hools 
B1_h(~mask)=B1match(~mask);
%save the B1 map in SPGR space
dtiWriteNiftiWrapper(single(B1_h), xform, B1file);
B1=B1_h;

return

%%%old code
%%
% filter2=50./pixdim;
% 
% [f2] = makegaussian3d(filter2,[0.5 0.5 0.5],[0.25 0.25 0.25]);
% C2 = convn(tisuuemask,f2,'same');
% CC2=convn(tt,f2,'same');
% 
% tisuuemask2=C2>max(CC2(:))*.3 & C1<max(CC1(:)).*0.6;
% [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask2));
% w2 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter2,[]);
% tmp2(find(tisuuemask2))=w2;
% 
% %%
% filter3=60./pixdim;
% 
% [f3] = makegaussian3d(filter3,[0.5 0.5 0.5],[0.25 0.25 0.25]);
% C3 = convn(tisuuemask,f3,'same');
% CC3=convn(tt,f3,'same');
% 
% tisuuemask3=C3>max(CC3(:))*.3 & C2<max(CC2(:)).*0.4;
% [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask3));
% w3 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter3,[]);
% tmp3(find(tisuuemask3))=w3;
% 
% 
% 
% %%
% filter4=70./pixdim;
% 
% [f4] = makegaussian3d(filter4,[0.5 0.5 0.5],[0.25 0.25 0.25]);
% C4 = convn(tisuuemask,f4,'same');
% CC4=convn(tt,f4,'same');
% 
% tisuuemask4=C4>max(CC4(:))*.3 & C3<max(CC3(:)).*0.4;
% [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask4));
% w4 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter4,[]);
% tmp4(find(tisuuemask4))=w4;
% 
% 
% 
% %%
% filter5=120./pixdim;
% 
% [f5] = makegaussian3d(filter5,[0.5 0.5 0.5],[0.25 0.25 0.25]);
% C5 = convn(tisuuemask,f5,'same');
% CC5=convn(tt,f5,'same');
% 
% tisuuemask5=C5>max(CC5(:))*.3 & C4<max(CC4(:)).*0.4;
% [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask5));
% w5 = localregression3d(x,y,z,B1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter4,[]);
% tmp5(find(tisuuemask5))=w5;
% 
% 
% %%
% tt4=(tmp1+tmp2+tmp3+tmp4+tmp5)./(tisuuemask1+tisuuemask2+tisuuemask3+tisuuemask4+tisuuemask5);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% area0=area;
% filterS=FS
% step=2;
% while i<intep 
%     
%      tisuuemask0=(logical(tmp1));
%     if length(find(tisuuemask0==0))==0 %this we happen when we are dwon    
%        %then we are dwon
%         i=intep+1;
%     end
%     
%     i=i+1 %other wise wwe will incrise the rigion we used to fit
%     filterS=filterS+step;
%      filter1=(filterS)./pixdim;
%         Farea1=filter1(1)*filter1(2)*filter1(3);
%         [f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);
%         CC1=convn(tt,f1,'same');
%         area=area0*Farea/Farea1
%     
%       if area<0.2 && filterS<60%no point to estimate if you don't have at list 20% of
%             % the data
%             area=0.2;
%             step=6;
%         end
%    if area<0.1 && filterS>60 %no point to estimate if you don't have at list 20% of
%             % the data but if we got so far and the fillter is big lets
%             % just finish and get it over with
%             area=0.1;
%             step=6;
%         end
%     C11 = convn(tisuuemask,f1,'same');
%     
%     tisuuemask1=C11>(max(CC1(:)).*area) & ~tisuuemask0;
%     if  length(find(tisuuemask1~=0))>0
%         [x y z]=ind2sub(size(tisuuemask),find(tisuuemask));
%         [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));
%         wn = localregression3d(x,y,z,tmp1(find(tisuuemask)),(x0),(y0),(z0),[],[],filter1,[]);
%         tmp1(find(tisuuemask1))=wn;
%     
%         end
% 
% end
% 




%% let try to interpolate and fit the B1 with the original B1 smooth model and the mesuraments we have
% we will alow some change in the filter size so we will get more covarage
% tic;
% 
% filter1=FS./pixdim;
% Farea=filter1(1)*filter1(2)*filter1(3);
% %add=tisuuemask & ~tisuuemask1;
% %tmp1(add)=B1(add);
% tisuuemask2=tisuuemask1;
% filterS=FS;
%          i=0;
% while i<intep
%     i=i+1
%     tisuuemask0=(logical(tmp1));
%     if length(find(tisuuemask0==0))==0
%       
%                 %then we are dwon
%         i=intep+1;
%     end
%     C11 = convn(tisuuemask0,f1,'same');
% 
%     tisuuemask1=C11>(max(CC1(:)).*area) & ~tisuuemask2;
%     if length(find(tisuuemask1~=0))>0
%         tisuuemask2=tisuuemask1+tisuuemask2;
%         [x y z]=ind2sub(size(tisuuemask),find(tisuuemask0));
%         [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));
%         wn = localregression3d(x,y,z,tmp1(find(tisuuemask0)),(x0),(y0),(z0),[],[],filter1,[]);
%         tmp1(find(tisuuemask1))=wn;
%     else
% 
%         filterS=filterS+2
%         filter1=(filterS)./pixdim;
%         %Farea1=filter1(1)*filter1(2)*filter1(3);
%         [f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);
%         CC1=convn(tt,f1,'same');
%         %area=area*Farea/Farea1
%         %if area<0.2 %no point to estimate if you don't have at list 20% of
%         %the data
%         %area=0.2;
%         %end
% 
%         if filterS>FS.*1.25
%             i=intep+1;
%         end
% 
%     end
% end
% toc
% %% let try to interpolate and fit the B1 with the original B1 smooth model and the mesuraments we have
% % we will alow some change in the filter size and pracent area that has data so we will get more covarage
% filter1=FS./pixdim;
% Farea=filter1(1)*filter1(2)*filter1(3);
% %add=tisuuemask & ~tisuuemask1;
% %tmp1(add)=B1(add);
% %tisuuemask2=tisuuemask1;
% filterS=FS;
% step=6;
%  area0=area;     i=0;
% tic;
% while i<intep
%     i=i+1
%     tisuuemask0=(logical(tmp1));
%     if length(find(tisuuemask0==0))==0
%         
%        %then we are dwon
%         i=intep+1;
%     end
%     C11 = convn(tisuuemask0,f1,'same');
%     
%     tisuuemask1=C11>(max(CC1(:)).*area) & ~tisuuemask2;
%     if  length(find(tisuuemask1~=0))>0
%         tisuuemask2=tisuuemask1+tisuuemask2;
%         [x y z]=ind2sub(size(tisuuemask),find(tisuuemask0));
%         [x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));
%         wn = localregression3d(x,y,z,tmp1(find(tisuuemask0)),(x0),(y0),(z0),[],[],filter1,[]);
%         tmp1(find(tisuuemask1))=wn;
%     else
%         filterS=filterS+step
%         filter1=(filterS)./pixdim;
%         Farea1=filter1(1)*filter1(2)*filter1(3);
%         [f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);
%         CC1=convn(tt,f1,'same');
%         area=area0*Farea/Farea1
%           if area<0.2 && filterS<60%no point to estimate if you don't have at list 20% of
%             % the data
%             area=0.2;
%             step=6;
%         end
%         if area<0.1 && filterS>60 %no point to estimate if you don't have at list 20% of
%             % the data but if we got so far and the fillter is big lets
%             % just estimate and get it over with
%             area=0.1;
%             step=6;
%         end
%         
%         
%       
%         
%         if filterS>200
%             i=intep+1;
%         end
%         
%     end
% end
% toc
% 
% save('tmp1B1map','tmp1');
% 
% %%
%     tisuuemask0=(logical(tmp1));
% if length(find(tisuuemask0==0))==0
%     
% 
%     end