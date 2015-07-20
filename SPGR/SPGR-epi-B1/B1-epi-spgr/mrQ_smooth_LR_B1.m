function mrQ=mrQ_smooth_LR_B1(mrQ,smoothnessVal)
%% load  the fit information


% load the fit parameeters
load  (mrQ.B1.logname)

B1=readFileNifti(mrQ.B1.epiFitFileName);
pixdim=B1.pixdim;   xform=B1.qto_xyz;
B1=double(B1.data);

resnormMap=readFileNifti(mrQ.B1.resnormFileName);   resnormMap=double(resnormMap.data);
UseVoxNMap=readFileNifti(mrQ.B1.NvoxFileName);     UseVoxNMap=double(UseVoxNMap.data);


%% we will smooth and interpulate/exstapulate the values to have solution in every location
%the fit error resnorm (sum of error) devided by number of voxels
errMap=(resnormMap./UseVoxNMap);

tissuemask=resnormMap>0;
% don't use the misfit location
tissuemask=tissuemask & errMap< prctile(errMap(tissuemask),95);



%% II fit local regressions
% we will smooth the B1 map by local regresiions
%the fit is done in tree steps
%1.we first estimate the voxel that have can be estimate with great confidance
%(<45% ) of the area under the filter
%local information

%2. we fit  the others that are possible but with less confidence with
%more smooth (bigger filter)

%3. the one that are out of reach for local regration will be fitted by
%global polynomyial along the full B1 space (like a fillter along the all B1 map)

B1Fit=zeros(size(tissuemask));


sz=size(tissuemask);
tt=ones(sz);

%%%
%1. we check the local information by comparing to covariance of
%the filter with a all ones image (full information).

area=0.45;
%filter size
FS=30;
filter1=FS./pixdim;
[f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);

%define the available coverage
C1 = convn(tissuemask,f1,'same');

%define the maximal caverage
CC1=convn(tt,f1,'same');

%the voxel that we will use
tissuemask1=C1>max(CC1(:)).*area;

%where there are B1 estimation (x,y,z location)
[x y z]=ind2sub(size(tissuemask),find(tissuemask));

%where we will find the smooth B1  estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tissuemask),find(tissuemask1));

%local regression
w1 = localregression3d(x,y,z,B1(find(tissuemask)),(x0),(y0),(z0),[],[],filter1,[]);

%save the result
B1Fit(find(tissuemask1))=w1;




%% exstrapulate by smooth surfaces

if notDefined('smoothnessVal')
    smoothnessVal=5; % this is a value for gridfit see inside.
end


B1Fit_S=B1Fit;

%loop over Z slices

[XI YI]=meshgrid(1:size(B1Fit,1),1:size(B1Fit,2));

for  jj=1:size(B1Fit,3)
    
    tmp=B1Fit(:,:,jj);
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.2 && length(find(tmp>0))>1000);
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth version of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        B1Fit_S(:,:,jj)=ZI;
        
        clear ZI
    end;
    
end;

%%

B1Fit_S(B1Fit_S<0)=0;
%%  %loop over  x slices


[XI YI]=meshgrid(1:size(B1Fit_S,2),1:size(B1Fit_S,3));

for  jj=1:size(B1Fit_S,1)
    
    tmp=squeeze(B1Fit_S(jj,:,:));
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
        
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth version of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
       
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
       
        B1Fit_S(jj,:,:)=ZI;

        
        
        clear ZI
    end;
    
end;


B1Fit_S(B1Fit_S<0)=0;

%%  %loop over  y slices

[XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,3));

for  jj=1:size(B1Fit_S,2)
    
    tmp=squeeze(B1Fit_S(:,jj,:));
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
        
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth version of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        B1Fit_S(:,jj,:)=ZI;
        
        clear ZI
    end;
    
end;
B1Fit_S(B1Fit_S<0)=0;


%% calculate if the smoothing intruduce a constant bias, and correct for it.
%remove values that would cause correction value to be 0 or inf.
tissuemask(B1Fit_S==0)=0;
tissuemask(isnan(B1(:)./B1Fit_S(:)))=0;
tissuemask(isinf(B1(:)./B1Fit_S(:)))=0;

Cal=median(B1(tissuemask)./B1Fit_S(tissuemask));
 B1Fit_S=B1Fit_S.*Cal;



%% SAVE the result smooth B1 map
outDir = mrQ.spgr_initDir; % check that this is right

B1epiFileName=fullfile(outDir,['B1epi_map.nii.gz']);
dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1epiFileName);

mrQ.B1.epiFileName=B1epiFileName;
save(mrQ.name,'mrQ');

%% clear the fit tmp files

%eval(['! rm ' opt.name '*']);