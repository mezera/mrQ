function mrQ=mrQ_smooth_LR_B1(mrQ,outDir,smoothnessVal)
% mrQ=mrQ_smooth_LR_B1(mrQ,smoothnessVal)
%
% This function smoothes the B1 fit in EPI space, before returning to SPGR
% space. It uses the selected value for smoothing, as it loops over the
% slices in each of the X, Y and Z directions.
%
% ~INPUTS~
%              mrQ: The mrQ structure
%    smoothnessVal: A scalar value for the smoothness of the grid. Default
%                          is 5, and was selected through trial and error.  
%                          See gridfit.m for additional information.
%
% ~OUTPUTS~
%              mrQ: The updated mrQ structure.
%
% See also: mrQ_build_epi2SPGR_B1.m and gridfit.m.
%
% Edited by A.M 2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%  2015
%
%

%% I. Load the fit information

% load the fit parameters
load(mrQ.B1.logname)

B1=readFileNifti(mrQ.B1.epiFitFileName);
pixdim=B1.pixdim;   xform=B1.qto_xyz;
B1=double(B1.data);

resnormMap=readFileNifti(mrQ.B1.resnormFileName);   resnormMap=double(resnormMap.data);
UseVoxNMap=readFileNifti(mrQ.B1.NvoxFileName);     UseVoxNMap=double(UseVoxNMap.data);

% We will smoothe and interpolate/extrapolate the values to have a solution 
% at every location.

%the fit error resnorm (sum of error) divided by number of voxels
errMap=(resnormMap./UseVoxNMap);

tissuemask=resnormMap>0;
% don't use the misfit location
tissuemask=tissuemask & errMap< prctile(errMap(tissuemask),95);
%tissuemask=tissuemask & errMap<0.15;
%% II. Fit local regressions
% We will smoothe the B1 map by using local regressions.

% The fit is done in three steps:
%    1. We first estimate the voxels that can be estimated with great confidence
%       (<45% of the area under the filter)
%    2. We fit the other voxels that are possible (but with less confidence) 
%       with a smoother (bigger) filter
%    3. The voxels that are out of reach for local regression will be fitted 
%       by global polynomial along the full B1 space (like a filter along 
%       the entire B1 map)

B1Fit=zeros(size(tissuemask));

sz=size(tissuemask);
tt=ones(sz);

%
% 1. We check the local information by comparing the covariance of
%   the filter that of an all-ones image (full information).

area=0.45; 

%filter size
FS=30;

filter1=FS./pixdim;

[f1] = makegaussian3d(filter1,[0.5 0.5 0.5],[0.25 0.25 0.25]);

%Define the available coverage
C1 = convn(tissuemask,f1,'same');

%Define the maximal coverage
CC1=convn(tt,f1,'same');

%The voxels that we will use:
tissuemask1=C1>max(CC1(:)).*area;

%Where there is B1 estimation (x,y,z location)
[x y z]=ind2sub(size(tissuemask),find(tissuemask));

%Where we will find the smooth B1 estimation (x,y,z location)
[x0 y0 z0]=ind2sub(size(tissuemask),find(tissuemask1));

%Local regression
w1 = localregression3d(x,y,z,B1(find(tissuemask)),(x0),(y0),(z0),[],[],filter1,[]);

%Save the result
B1Fit(find(tissuemask1))=w1;

%% III. Extrapolate by smooth surfaces

if notDefined('smoothnessVal')
    smoothnessVal=5; % this is a value for gridfit; see inside.
end


B1Fit_S=B1Fit;

%% IIIa. Loop over Z slices

[XI YI]=meshgrid(1:size(B1Fit,1),1:size(B1Fit,2));

for  jj=1:size(B1Fit,3)
    
    tmp=B1Fit(:,:,jj);
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.2 && length(find(tmp>0))>1000);
            %          This if makes aure we only extrapolate on data
            %          that's large enough within the scan. but it's not
            %          done beautifully, and maybe a better way would be to
            %          check the ratio of B1 voxels in alice and the number of
            %          voxels in a brainmmask in that same slice (ant not
            %          the number of voxel in that slice regardless of
            %          bain, as it is done now).
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
        if  ~isempty(find(isnan(ZI),1)) % we might get NaNs in the edges
            ZIt = griddata(xg,yg,zg,XI,YI,'v4');
            ZI(isnan(ZI))=ZIt(isnan(ZI));
        end
                        ZI(isnan(ZI))=0;

        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        B1Fit_S(:,:,jj)=ZI;
        
        clear ZI
    end;
    
end;

B1Fit_S(B1Fit_S<0)=0;

%%  IIIb. Loop over X slices

[XI YI]=meshgrid(1:size(B1Fit_S,2),1:size(B1Fit_S,3));

for  jj=1:size(B1Fit_S,1)
    
    tmp=squeeze(B1Fit_S(jj,:,:));
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
        
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
        if  ~isempty(find(isnan(ZI),1))% we might get NaNs in the edges
            ZIt = griddata(xg,yg,zg,XI,YI,'v4');
            ZI(isnan(ZI))=ZIt(isnan(ZI));
        end
        
               ZI(isnan(ZI))=0;

        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
       
        B1Fit_S(jj,:,:)=ZI;
        
        clear ZI
    end;
    
end;

B1Fit_S(B1Fit_S<0)=0;

%%  IIIc. Loop over Y slices

[XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,3));

for  jj=1:size(B1Fit_S,2)
    
    tmp=squeeze(B1Fit_S(:,jj,:));
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
        
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
        ZI = griddata(xg,yg,zg,XI,YI);
        if  ~isempty(find(isnan(ZI),1)) % we might get NaNs in the edges
            ZIt = griddata(xg,yg,zg,XI,YI,'v4');
            ZI(isnan(ZI))=ZIt(isnan(ZI));
        end
        
        ZI(isnan(ZI))=0;
        
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        B1Fit_S(:,jj,:)=ZI;
        
        clear ZI
    end;
    
end;

B1Fit_S(B1Fit_S<0)=0;

%% IV. Calculate if the smoothing introduces a constant bias, and correct for it.
%remove values that would cause correction value to be 0 or inf.
tissuemask(B1Fit_S==0)=0;
tissuemask(isnan(B1(:)./B1Fit_S(:)))=0;
tissuemask(isinf(B1(:)./B1Fit_S(:)))=0;

Cal=median(B1(tissuemask)./B1Fit_S(tissuemask));
 B1Fit_S=B1Fit_S.*Cal;

%% V. SAVE the resulting smooth B1 map

B1epiFileName=fullfile(outDir,['B1epi_map.nii.gz']);
dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1epiFileName);

mrQ.B1.epiFileName=B1epiFileName;
save(mrQ.name,'mrQ');

%% clear the fit tmp files

%eval(['! rm ' opt.name '*']);