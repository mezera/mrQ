function mrQ=mrQ_build_epi2SPGR_B1(mrQ,OutDir,TargetT1file,B1FileName,smoothnessVal)
% function mrQ=mrQ_build_epi2SPGR_B1(mrQ,B1FileName,smoothnessVal)
%
% In this function, we will smoothe and interpolate/extrapolate the values
% to have a solution at every location. The EPI image is warped into SPGR
% space, and then a similar smoothing mechanism is employed (using the
% selected value for smoothing, then looping over the slices in each of the
% X, Y and Z directions).
%
% ~INPUTS~
%              mrQ: The mrQ structure
%       B1FileName: Specify where the B1 inhomogeneity map exists. If you 
%                      leave it empty, it will use the B1 file from the 
%                      data directory.
%    smoothnessVal: Value for the smoothness of the grid. Default is 5.
%
% ~OUTPUTS~
%              mrQ: The updated mrQ structure
%
% SEE ALSO: mrQ_smooth_LR_B1.m
%
% Edited A.M 2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%

%% I. Load files and set parameters


if ~isfield(mrQ.Ants_Info,'B1_epi_spgr')
     mrQ.Ants_Info=mrQ_ANTS_warp_EPI2SPGR(mrQ.Ants_Info,TargetT1file,OutDir,mrQ.B1.epiFileName);
end

T1=readFileNifti(mrQ.Ants_Info.T1_epi_spgr);
B1=readFileNifti(mrQ.Ants_Info.B1_epi_spgr);

xform=T1.qto_xyz;
mask=logical(T1.data);
B1Fit_S=zeros(size(mask));

mask(B1.data==0)=0;
mask(isinf(B1.data))=0;
mask(isnan(B1.data))=0;

B1Fit_S(mask)=B1.data(mask);
% clear T1 B1 mask;

clear T1;

if notDefined('smoothnessVal')
smoothnessVal=5; % this is a value for gridfit; see inside.
end

%% IIa. Loop over Z slices

[XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,2));

    for  jj=1:size(B1Fit_S,3)
        
        tmp=B1Fit_S(:,:,jj);
        
        % Check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            
            % Estimate a smooth version of the data in the slice
            %
            % (For original code, see: 
            %     Noterdaeme et al.
            %     "Intensity correction with a pair of spoiled gradient
            %       recalled echo images".
            %     Phys. Med. Biol. 54 3473-3489 (2009)
            %                                         )
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
            ZI = griddata(xg,yg,zg,XI,YI);
            
                if  ~isempty(find(isnan(ZI),1)) % We might get NaNs in the edges
                    ZIt = griddata(xg,yg,zg,XI,YI,'v4');
                    ZI(isnan(ZI))=ZIt(isnan(ZI));
                end
                
            % Put the resulting gain in the 3D gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            B1Fit_S(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
    B1Fit_S(B1Fit_S<0)=0;

%% IIb. Loop over X slices

    [XI YI]=meshgrid(1:size(B1Fit_S,2),1:size(B1Fit_S,3));

    for  jj=1:size(B1Fit_S,1)
        
        tmp=squeeze(B1Fit_S(jj,:,:));
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth version of the data in the slice
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
            ZI = griddata(xg,yg,zg,XI,YI);
            if ~isempty(find(isnan(ZI),1)) % we might get nan in the edges
                ZIt = griddata(xg,yg,zg,XI,YI,'v4');
                ZI(isnan(ZI))=ZIt(isnan(ZI));
            end
            % Put the resulting gain in the 3D gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            B1Fit_S(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
        B1Fit_S(B1Fit_S<0)=0;

    
%% IIc. Loop over Y slices
    
    [XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,3));

    for  jj=1:size(B1Fit_S,2)
        
        tmp=squeeze(B1Fit_S(:,jj,:));
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth version of the data in the slice
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
            ZI = griddata(xg,yg,zg,XI,YI);
            
                if  ~isempty(find(isnan(ZI),1)) % we might get NaNs in the edges
                    ZIt = griddata(xg,yg,zg,XI,YI,'v4');
                    ZI(isnan(ZI))=ZIt(isnan(ZI));
                end
                
            % Put the resulting gain in the 3D gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            B1Fit_S(:,jj,:)=ZI;
            
            clear ZI
        end;
        
    end;   
        B1Fit_S(B1Fit_S<=0)=0;

%% III. Calculate if the smoothing introduces a constant bias, and correct for it.

% Original bias clearing (from mrQ_smooth_LR_B1.m)
mask(B1Fit_S==0)=0;

mask(  isnan( B1.data(:)./B1Fit_S(:) )  )=0;
mask(  isinf( B1.data(:)./B1Fit_S(:) )  )=0;

Cal=median(B1.data(mask)./B1Fit_S(mask));
B1Fit_S=B1Fit_S.*Cal; 

% B1Fit_S(B1Fit_S<=0)=1;
% B1Fit_S(isnan(B1Fit_S))=1;
% B1Fit_S(isinf(B1Fit_S))=1;

B1Fit_S(B1Fit_S<=0)=eps;
B1Fit_S(isnan(B1Fit_S))=eps;
B1Fit_S(isinf(B1Fit_S))=eps;


%% IV. Save
%outDir = mrQ.spgr_initDir; 

if notDefined('B1FileName')
        B1FileName=fullfile(OutDir,'B1_Map.nii.gz');
end

        dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1FileName)
        
        mrQ.B1FileName=B1FileName;
      save(mrQ.name,'mrQ');
      
  fprintf(['Done fitting B1 map. B1 file is saved: '   mrQ.B1FileName        '  \n']);

