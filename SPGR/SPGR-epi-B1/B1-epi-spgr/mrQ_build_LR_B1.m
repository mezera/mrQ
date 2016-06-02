function B1files=mrQ_build_LR_B1(B1logname,outDir)
% function mrQ=mrQ_build_LR_B1(mrQ)
%
%  In this function, the estimates of each single voxel are joined
%together to form a 3D B1 image. It uses the data generated for the single
%voxels from a previous function in the B1_LR pipeline (mrQ_B1FitMask.m),
%as stored in the mrQ structure.
%
% Accepts the mrQ structure as an INPUT and returns an updated mrQ
% structure as the OUTPUT.
%
%
% See also: mrQ_B1FitMask, mrQ_PD_LRB1SPGR_GridParams
%
% Edited by A.M 2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%

%% I. Load the fit parameters
        load(B1logname)
B1files.logname=B1logname;
    %% II. Find where the fitted values are in 3D space

        BM=readFileNifti(opt.tisuuemaskFile);
        xform=BM.qto_xyz;
        BM=logical(BM.data);
        BM= true(size(BM));

     loc=find(BM);
%% III. Loop over voxels, and load the fitted values

jumpindex=opt.jumpindex; % number of voxels in each file
jobindexs=1:ceil(opt.N_Vox2Fit/jumpindex);  % the numbers of saved file each with jumpindex 

B1Fit=zeros(size(BM));
resnormMap=B1Fit;
UseVoxNMap=B1Fit;

for ii=1:length(jobindexs);
        % the number of the box in the file (also part of the saved file name)
    
    st=1 +(jobindexs(ii)-1)*jumpindex;
    ed=st+jumpindex-1;
    
    %check that this box has brain data
    if ed>opt.N_Vox2Fit, ed=opt.N_Vox2Fit;end;
    FileName=[ opt.name '_' num2str(st) '_' num2str(ed)];
    load(FileName);
    FitVoxel=loc(st:ed);
    B1Fit(FitVoxel(skip~=1))=B1((skip~=1));
    UseVoxNMap(FitVoxel(skip~=1))=UseVoxN((skip~=1));
    resnormMap(FitVoxel(skip~=1))=resnorm((skip~=1));
end
   
%% IV. Save 3D results as NIfTI
%outDir = mrQ.spgr_initDir; %check that this is right

% the B1 map
 B1epiFitFileName=fullfile(outDir,['B1_LR.nii.gz']);
dtiWriteNiftiWrapper(single(B1Fit),xform,B1epiFitFileName);

% the fitting error
B1resnormFileName=fullfile(outDir,['B1resnorm_LR.nii.gz']);
dtiWriteNiftiWrapper(single(resnormMap),xform,B1resnormFileName);

% The number of estimations for each voxel. 
% Note that each time we have fit several voxels, and that there is overlap 
% between those fits.

B1NvoxFileName=fullfile(outDir,['B1_Nvox_LR.nii.gz']);
dtiWriteNiftiWrapper(single(UseVoxNMap),xform,B1NvoxFileName);

%  Record the path to the B1 fit
B1files.epiFitFileName=B1epiFitFileName;
B1files.resnormFileName=B1resnormFileName;
B1files.NvoxFileName=B1NvoxFileName;


%%
% We will smoothe and interpolate/extrapolate the values to have solutions
% in every location. The fit error rensorm (sum of error) divided by the number
% of voxels errMap=(resnormMap./UseVoxNMap);
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
%             % Noterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
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
%  B1epiFileName=fullfile(AnalysisInfo.outDir,['B1epi_map.nii.gz']);
% dtiWriteNiftiWrapper(single(B1Fit_S),xform,B1epiFileName);
% 
% AnalysisInfo.B1epiFileName=B1epiFileName;
%         save(AnalysisInfo.name,'AnalysisInfo');

%% clear the fit tmp files

%clear the fit files
%eval(['! rm ' opt.name '*']);