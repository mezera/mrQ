function  [Gain1]=mrQ_smmothL_GainPhantoms(T1,M0,outDir,xform,mrQ,Gainfile,SPGRResidfile,mask);

%
%
% Perform smoth the a coil gaincorection to phantom for smooth space
%
% INPUTS:
%   
% OUTPUTS
% Gain -smooth B1 in SPGR space
%AnalysisInfo -info structure updated after the registration



%we smoth the SEIR B1 map by local and then feel the gap by global
%regrigions and then register back the smooth B1 to the SPGR space by Ants
%softwere

% See Also:
% mrQfit_T1M0_ver2

 
if (~exist('degree','var')|| isempty(degree)),
    degree=3;
end;

if notDefined('SPGRResidfile')
    SPGRResidfile=fullfile(outDir,['lsqT1PDresnorm_last.nii.gz']);
end
%%  B1 mask
% mask the B1 value that we don't like to use 
%1. more then 50% T1 variation between SEIR and SPGR is probabaly miss registration
%2. B1 in the edge of the serech space is point that didn't converge in the B1
%fit
%3.voxel with big fit residual for B1 or SEIR T1 are also out

%load the B1 fit

if notDefined('mask')
    
       mask = mrQ.hoginiues_mask;   
             
end

 tisuuemask =readFileNifti(mask);
 
tisuuemask=double(tisuuemask.data);
tisuuemask_forcoils=zeros(size(tisuuemask));
             tisuuemask_forcoils(tisuuemask==1)=1;
tisuuemask_forcoils=logical(tisuuemask_forcoils);
%load the data that was us to the B1 fit


%load the B1 fit residual
fitResid=readFileNifti(SPGRResidfile);
pixdim=fitResid.pixdim;

fitResid= fitResid.data;
%load the T1 fit residuals



%no fot to big or small B1 (50%)
tisuuemask2=  tisuuemask_forcoils & fitResid<prctile(fitResid(find(fitResid)),97) ;


% no for big residual

   
tisuuemask2=tisuuemask2 & T1>prctile(T1(tisuuemask2),1) & T1<prctile(T1(tisuuemask2),99);

tisuuemask2=tisuuemask2 & M0>prctile(M0(tisuuemask2),1) & M0<prctile(M0(tisuuemask2),99);

%%
Gain=zeros(size(M0));
Gain(tisuuemask2)=M0(tisuuemask2);

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
C1 = convn(tisuuemask2,f1,'same');

%define the maxsimal caverage
CC1=convn(tt,f1,'same');

%the voxel that we will use
tisuuemask1=C1>max(CC1(:)).*area;
%tisuuemask1=tisuuemask1 & tisuuemask;
%where there are gain estimation (x,y,z location)
[x y z]=ind2sub(size(tisuuemask),find(tisuuemask2));

%were we will find the smooth gain  estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));

%local regrision
w1 = localregression3d(x,y,z,Gain(find(tisuuemask2)),(x0),(y0),(z0),[],[],filter1,[]);

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
C1 = convn(tisuuemask2,f1,'same');
%where there are B1 estimation

%define the maxsimal caverage
CC1=convn(tt,f1,'same');

%the voxel that we will use
tisuuemask1=C1>max(CC1(:)).*area;
%tisuuemask1=tisuuemask1 & tisuuemask;

%where there are B1 estimation (x,y,z location)
[x y z]=ind2sub(size(tisuuemask),find(tisuuemask2));

%were we will find the smooth B1 estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tisuuemask),find(tisuuemask1));

%local regrision
w2 = localregression3d(x,y,z,Gain(find(tisuuemask2)),(x0),(y0),(z0),[],[],filter1,[]);

%save the result
tmp1(find(tisuuemask1))=w2;




%%% 3.   we can do global polynomial if we can't fit localy some locations    
Imsz1=size(tmp1);

%make the 3D polynials
[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),degree);

%fit the polynials coefitents to the smooth B1 map
[params,gains,rs] = fit3dpolynomialmodel(tmp1,(tmp1>0),degree);

%reshape from a vector to 3D map
B1match = reshape(Poly1*params(:),Imsz1);

%find where the B1 value are missing
mask1=logical(tmp1>0);

%fill the holls
tmp1(~mask1)=B1match(~mask1);


%%% last check
%we might get  scaling effect when the smoothing is not perefect , scale it
%back so the mean of the smooth map and the original is the same




%%% save the smooth B1 map
dtiWriteNiftiWrapper(single(tmp1), xform, Gainfile);



Gain1=tmp1;

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