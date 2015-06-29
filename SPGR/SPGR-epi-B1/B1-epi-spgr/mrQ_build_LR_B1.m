function mrQ=mrQ_build_LR_B1(mrQ)
% join the estimate of each single voxel to a 3D B1 image.

%%

% load the fit parameeters
    load  (mrQ.B1.logname)
      



%% find where are the fitted value in 3D space

    BM=readFileNifti(opt.tisuuemaskFile);
    xform=BM.qto_xyz;
    BM=logical(BM.data);
    BM=logical(ones(size(BM)));

 loc=find(BM);
%% loop over voxels load the fited values

jumpindex=opt.jumpindex; % number of voxel in each file
jobindexs=1:ceil(opt.N_Vox2Fit/jumpindex);  % the numbers of saved file each with jumpindex 

B1Fit=zeros(size(BM));
resnormMap=B1Fit;
UseVoxNMap=B1Fit;
for ii=1:length(jobindexs);
    
    % the number of the box in the file (also part of the saved file name)
    
    st=1 +(jobindexs(ii)-1)*jumpindex;
    ed=st+jumpindex-1;
    
    %cheack that this box have brain data
    if ed>opt.N_Vox2Fit, ed=opt.N_Vox2Fit;end;
    FileName=[ opt.name '_' num2str(st) '_' num2str(ed)];
    load(FileName);
    FitVoxel=loc(st:ed);
    B1Fit(FitVoxel(skip~=1))=B1((skip~=1));
    UseVoxNMap(FitVoxel(skip~=1))=UseVoxN((skip~=1));
    resnormMap(FitVoxel(skip~=1))=resnorm((skip~=1));
end
   
%% save 3D results as nifti
outDir = mrQ.spgr_initDir; % check that this is right

% the B1 map
 B1epiFitFileName=fullfile(outDir,['B1_LR.nii.gz']);
dtiWriteNiftiWrapper(single(B1Fit),xform,B1epiFitFileName);
% the fiting error
B1resnormFileName=fullfile(outDir,['B1resnorm_LR.nii.gz']);
dtiWriteNiftiWrapper(single(resnormMap),xform,B1resnormFileName);
% the number of estimation for each voxel. not that we have fit each time
% several voxel and there were overlap betwwn those fits.
B1NvoxFileName=fullfile(outDir,['B1_Nvox_LR.nii.gz']);
dtiWriteNiftiWrapper(single(UseVoxNMap),xform,B1NvoxFileName);

%  record the path to the B1 fit
mrQ.B1.epiFitFileName=B1epiFitFileName;
mrQ.B1.resnormFileName=B1resnormFileName;
mrQ.B1.NvoxFileName=B1NvoxFileName;

      save(mrQ.name,'mrQ');

    
    
    
    
    
    
    
    
%%
% %% we will smooth and interpulate/exstapulate the values to have solution in every location
% %the fit error rensorm (sum of error) devided by number of voxels
% errMap=(resnormMap./UseVoxNMap);
% 
% mask=resnormMap>0;
% % don't use the misfit location
% mask=mask & errMap< prctile(errMap(mask),95);
% 
% 
% if notDefined('smoothnessVal')
% smoothnessVal=5; % this is a value for gridfit see inside.
% end
% %% Local regration with biger window
% 
% 
% %%
% B1Fit_S=zeros(size(B1Fit));
% 
% [XI YI]=meshgrid(1:size(B1Fit,1),1:size(B1Fit,2));
% 
%  %loop over slices
%     for  jj=1:size(B1Fit,3)
%         
%         tmp=B1Fit(:,:,jj);
%         
%         %check that there is data in the slice
%         wh=find(tmp>0);
%         if  length(find(tmp>0))>100;
%             %find location of data
%             [x,y] = ind2sub(size(tmp),wh);
%             z=double(tmp(wh));
%             % estimate a smooth vertion of the data in the slice for original code see:
%             % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
%             
%             [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',smoothnessVal);
%             ZI = griddata(xg,yg,zg,XI,YI);
%             % put the result gain in the 3D gain image and fix orientation
%             ZI=rot90(ZI);
%             ZI = flipdim(ZI,1);
%             B1Fit_S(:,:,jj)=ZI;
%             
%             clear ZI
%         end;
%         
%     end;
% 
% 
% 
%  B1epiFileName=fullfile(AnalysisInfo.outDir,['B1epi_map.nii.gz']);
% dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1epiFileName);
% 
% AnalysisInfo.B1epiFileName=B1epiFileName;
%         save(AnalysisInfo.name,'AnalysisInfo');

%% clear the fit tmp files


%clear the fit files
%eval(['! rm ' opt.name '*']);