function mrQ=mrQ_B1_LR(mrQ)
% this function control the B1 fit.
% The function have two parts. 
% 1. fitting each voxel.
%     The fit is not independed to each voxel. A local regresion (LR) approch is taken.
%     When each voxel is fitted is close seraund voxel ae also used for the fit.
% 
% 2. Join all the voxels to B1 map.
% 	after join the voxels the fits we check for smoothnesswe  and exstrapulate to voxel that are missing in EPI space and then in SPGR space.

%% get ready to fit B1 and call to fit
if ~isfield(mrQ,'B1fit_done');
    mrQ.B1fit_done=0;
end

if ( mrQ.B1fit_done==0)
    
    
    % we build a mask for the voxel we like to fit with in epi space
% [mrQ.maskepi_File] = mrQ_B1FitMask(mrQ.Ants_Info.dirAnts,mrQ.spgr2epiAlignFile,mrQ.SEIR_epi_fitFile,mrQ.spgr_initDir);
% shai : there is not DirAnts in the Ants_Info, looks like you meant te
% directory where the aligned and warped files are.. 
[mrQ.maskepi_File] = mrQ_B1FitMask(mrQ.spgr_initDir,mrQ.spgr2epiAlignFile,mrQ.SEIR_epi_fitFile,mrQ.spgr_initDir);

    
    %define the fit parameters
    mrQ.B1.logname=mrQ_PD_LRB1SPGR_GridParams(mrQ);
    
    
    % call to fit (with or without grid)
    mrQ_fitB1LR_Call(mrQ.B1.logname,mrQ.SunGrid);
    
    
    %check that the fit is done (for SGE )  before you move on (while loop)
    mrQ. B1fit_done=mrQ_Gridcheack(mrQ.B1.logname,mrQ.SunGrid,3);
    
    save(mrQ.name,'mrQ');
    
end
    
    %%
    
    % build the grid B1 fits
    if isfield(mrQ,'B1Build_done');
    else
        mrQ.B1Build_done=0;
    end
    
    if (mrQ.B1Build_done==0 && mrQ.B1fit_done==1)
        fprintf('Build the B1 map form local fits           \n ');
        
        % join the estiations
        mrQ=mrQ_build_LR_B1(mrQ);
        
% smooth the B1 fit clear outlayer and exstarpulate for voxel that have no
% fits.
         mrQ=mrQ_smooth_LR_B1(mrQ);
         
         
        mrQ=mrQ_build_epi2SPGR_B1(mrQ);
                 mrQ.B1Build_done=1;
    save(mrQ.name,'mrQ');
        
    end    
            
 
    
    %% non linear Fit of  B1
    
%     %%% FIT B1 by lsq fit compare T1 SEIR(Res{1}.im) to the multi flip
%     % angle Res{3:end}.im % USE sge make the fit faster
%     
%     intM0= double(mean(M0(brainMask)));
%     flipAngles = [s(:).flipAngle];
%     tr = [s(:).TR];
%     
%     if(~all(tr == tr(1))), error('TR''s do not match!'); end
%     tr = tr(1);
%     
%     % Load or Create the tissuemask from all non-zero points
%     tisuuemaskFile=fullfile(outDir,['maskepiF.nii.gz']);
%     
%     if exist(tisuuemaskFile,'file') && ~clobber
%         tisuuemask_=readFileNifti(tisuuemaskFile);
%         tisuuemask_=logical(tisuuemask_.data);
%         
%     else
%         
%         tisuuemask_ =zeros(size(Res{1}.im));
%         
%         % Binarize the mask
%         tisuuemask_(find(Res{1}.im>10 & Res{3}.im>0))=1;
%         
%         % Create a logical array from the tissue mask to index the
%         % non-zero locations
%         tisuuemask_=logical(tisuuemask_);
%         
%         % Save the tissue mask
%         dtiWriteNiftiWrapper(single(tisuuemask_), SE_Xform, tisuuemaskFile);
%         
%     end;
%     
%     if clobber && (exist([outDir '/tmpSG'],'dir'))
%         % in the case we start over and there are old  fits we will
%         % deleat them
%         eval(['! rm -r ' outDir '/tmpSG']);
%     end
%     
%     % USE sge make the B1 fit faster
%     
%     
%     % This is lsq fit that uses the grid but you can make it not use
%     % SGE: see help inside mrQ_fitB1_LSQ
%     [B1 resNorm dd] = mrQ_fitB1_LSQ(Res, tisuuemask_, tr,flipAngles, outDir, intM0, SE_Xform, SunGrid, 1, [sub 'B1fit'],mrQ.proclus);
%     dtiWriteNiftiWrapper(single(B1), SE_Xform, B1epifile);
%     dtiWriteNiftiWrapper(single(resNorm), SE_Xform, B1epiResidfile);
%     
% %end
% 
% %% Move B1 from epi space to SPGR smooth and upsample
% 
% 
% 
% % Smooth the SEIR B1 map by local and then fill the gap by global
% % regretions and then register back the smoothed B1 to the SPGR space using ANTS
% [B1,AnalysisInfo]=mrQ_smmothL_B1(B1epifile,AlignFile,outDir,B1file,xform,[],AnalysisInfo,t1fileHM,SET1Fitfile,B1epiResidfile);
