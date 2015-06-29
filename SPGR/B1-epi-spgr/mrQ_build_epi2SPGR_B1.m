function mrQ=mrQ_build_epi2SPGR_B1(mrQ,B1FileName,smoothnessVal)

%%



if ~isfield(mrQ.Ants_Info,'RB_B1_epi_spgr')
     mrQ.Ants_Info=mrQ_RB_ANTS_warp_EPI2SPGR(mrQ.Ants_Info,mrQ.t1fileHM,mrQ.spgr_initDir,mrQ.B1.epiFileName);
%
end
T1=readFileNifti(mrQ.Ants_Info=.RB_T1_epi_spgr);
B1=readFileNifti(mrQ.Ants_Info=.RB_B1_epi_spgr);
xform=T1.qto_xyz;
mask=logical(T1.data);


B1Fit_S=zeros(size(mask));
B1Fit_S(mask)=B1.data(mask);
clear T1 B1 mask;



%%
%% we will smooth and interpulate/exctrapulate the values to have solution in every location


if notDefined('smoothnessVal')
smoothnessVal=5; % this is a value for gridfit see inside.
end


[XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,2));

 %loop over  zslices
    for  jj=1:size(B1Fit_S,3)
        
        tmp=B1Fit_S(:,:,jj);
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
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
    
    B1Fit_S(B1Fit_S<0)=0;
    
    %%
    
    [XI YI]=meshgrid(1:size(B1Fit_S,2),1:size(B1Fit_S,3));

 %loop over  x slices
    for  jj=1:size(B1Fit_S,1)
        
        tmp=squeeze(B1Fit_S(jj,:,:));
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
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

    
    %%
    
    [XI YI]=meshgrid(1:size(B1Fit_S,1),1:size(B1Fit_S,3));

 %loop over  y slices
    for  jj=1:size(B1Fit_S,2)
        
        tmp=squeeze(B1Fit_S(:,jj,:));
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmp(:))<1  && length(wh)>1000);
            
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
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
        B1Fit_S(B1Fit_S<=0)=0;

%% save
outDir = mrQ.spgr_initDir; % check that this is right

if notDefined(B1FileName)
        B1FileName=fullfile(AnalysisInfo.outDir,'B1_Map.nii.gz');
end

        dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1FileName)
        
        mrQ.B1FileName=B1FileName;
      save(mrQ.name,'mrQ');
      
              fprintf(['Done fitting B1 map. B1 file is saved: '   mrQ.B1FileName        '  \n');

