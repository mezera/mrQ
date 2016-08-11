function [AnalysisInfo]=mrQ_RB_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%
%[AnalysisInfo]=mrQ_RB_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%
% This function evaluates the ANTs affine linear registration from EPI to
% SPGR.  
%
% We apply the linear registration on the B1 map that was fitted in
% EPI space,  and bring it back to the SPGR space. The B1 will therefore be
% interpolated, but we are fine with that because we assume that B1 is a
% smooth parameter in the imaging space. The SPGR T1 is registered to the
% EPI warp version that was calculated before by
% mrQ_NLANTS_warp_SPGR2EPI_RB.m. This ensures that the registration is
% between two images with the same contrast. 
%
% The ANTs will register the EPI T1 to the SPGR T1 using unix functions.
% The code assumes ANTs and ITK code are part of the computer unix sell
% path (./bashrc). The registration happens in two steps: First, define the
% parameters; and second, apply it and warp the images. We warp both the T1
% and B1.
%
% The process is documented in the AnalysisInfo structure.
%
% INPUTS:
%
%      AnalysisInfo:  A structure that keeps records of the T1 fit files.
%
%          t1fileHM:  The path to SPGR T1 file. 
%                         (This is the uncorrected fast fit T1 map) 
%
%            outDir:  The path to where files are read from and written to
%
%            B1file:  The path to the B1 file in EPI space 
%                         (This file was saved by mrQ_smooth_LR_B1.m)
%
% OUTPUTS:
%
%      AnalysisInfo:  The updated structure that keeps records of the T1 
%                          fit files. The ANTs registration files' name and 
%                          date will be recorded in it.
%
% SEE ALSO: mrQ_T1M0_Fit.m and mrQ_smooth_LR_B1.m
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%


%% I. We first need to make matlab work with the bashrc that we have
%
% This was found to be unnecessary, but if ANTs doesn't seem to work -  
% uncomment this section to add it to the path, don't forget to un-comment 
% the path reset at the end of the script:
% 
% orig_path = getenv('LD_LIBRARY_PATH');
% 
% colon_idx = strfind(orig_path, ':');
% 
% matlab_path = [];
% other_path = [];
% 
% for i=1:(length(colon_idx)-1)
% 
%     this_path = orig_path(colon_idx(i):colon_idx(i+1));
%     % This is part of the matlab path:
%     if strfind(this_path, 'matlab')
%         matlab_path = [matlab_path, ':', this_path];
% 
%         % Should go before the matlab bit of the path (all the rest):
%     else
%         other_path = [other_path, ':', this_path];
% 
%     end
% end
% 
% new_path = [other_path, ':', matlab_path];
% 
% setenv('LD_LIBRARY_PATH', new_path)

%%
 %% II. Run ANTs non-linear parameter to register SPGR to EPI

%name the saved parameter file
AnalysisInfo.RB_WARP_EPI_SPGR=fullfile(outDir,'WARP_EPI_SPGR_RB');
%AnalysisInfo.RBMI_WARP_EPI_SPGR=fullfile(outDir,'WARP_EPI_SPGR_RBMI');


%if ~(exist([AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz'],'file') && exist([AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'],'file') )

    % ANTs: make the ANTs comand
    cmANTS = ['xterm -e ANTS 3 -m CC[' t1fileHM ',' AnalysisInfo.T1_spgr_epi_RB ',1,2]  -o ' AnalysisInfo.RB_WARP_EPI_SPGR '.nii.gz --rigid-affine true'];
    %cmANTS=['xterm -e ANTS 3 -m MI[' t1fileHM ',' AnalysisInfo.T1_spgr_epi ',1,32]   -o ' AnalysisInfo.RBMI_WARP_EPI_SPGR '.nii.gz --rigid-affine true']
    % Run the command in unix and get back status and results:
    [status, result] = system(cmANTS,'-echo');
    if status ~=0
        cmANTS=['ANTS 3 -m CC[' t1fileHM ',' AnalysisInfo.T1_spgr_epi_RB ',1,2]  -o ' AnalysisInfo.RB_WARP_EPI_SPGR '.nii.gz --rigid-affine true'];
        [status, result] = system(cmANTS,'-echo')
    end
    


%% III. Name and save the SPGR B1 map in SPGR space
% The T1 file
AnalysisInfo.RB_B1_epi_spgr=fullfile(outDir,'Warp_B1_EPI2SPGR_RB.nii.gz');
%cmWarp=['xterm -e WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi ' '  AnalysisInfo.B1_epi_spgr ' -R '  B1file  ' ' AnalysisInfo.WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.WARP_EPI_SPGR 'Affine.txt'];

    % ANTs: make the ANTs command
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' B1file ' '  AnalysisInfo.RB_B1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~=0
        cmWarp=['WarpImageMultiTransform  3 ' B1file ' '  AnalysisInfo.RB_B1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
        % Run the command in unix and get back status and results:
        [status, result] = system(cmWarp);
    end

% The B1 file
AnalysisInfo.RB_T1_epi_spgr=fullfile(outDir,'Warp_T1_EPI2SPGR_RB.nii.gz');

    % ANTs: make the ANTs command
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi_RB ' '  AnalysisInfo.RB_T1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~= 0
        cmWarp=['WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi_RB ' '  AnalysisInfo.RB_T1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
        [status, result] = system(cmWarp,'-echo');
    end


%% IV. Reset the ld_path:
  % Reset the ld_path (unnecessary, see explanation at top of function)
% setenv('LD_LIBRARY_PATH', orig_path)

%% V. Save the updated AnalysisInfo file
infofile=fullfile(outDir,'AnalysisInfo.mat');
save(infofile,'AnalysisInfo');
