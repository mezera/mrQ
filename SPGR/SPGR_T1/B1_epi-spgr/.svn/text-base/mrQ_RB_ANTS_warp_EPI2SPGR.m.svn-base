function [AnalysisInfo]=mrQ_RB_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%
%[AnalysisInfo]=mrQ_RB_ANTS_warp_EPI2SPGR(AnalysisInfo,t1fileHM,outDir,B1file)
%

% # this function evaluate the ANTS affine linear reistration from epi to
% spgr.   we applay it on the B1 map that was fitted in epi space and by
% that bring it back to the SPGR space. the B1 will be interpulate
% therefore but we are fine with that becouse we assume that B1 is a smooth
% parameter in the imaging space.
% the SPGR T1 is register to it epi warp version that was cacluate before by
% mrQ_NLANTS_warp_SPGR2EPI_RB.m. this promise that the registration are between to images with the same contrast. 
%the ANTS will register the  epi T1 tp the spgr T1 using unixs functions. the code esume ANTS and ITK code are part of the computer unix sell path (./bachrc).
%the registration happen in two steps first define parameters and secound
%apply it and wrap the images. we Warps both the T1 and B1.
% the procces is document in the  AnalysisInfo structure. 
%
% INPUTS:
%
%    AnalysisInfo -  a structure that keep records of the T1 fit files
%
%   
%   t1fileHM  - the path to SPGR T1 file (this is uncorected fast fit T1 map). 
%
%
%   outDir - the path to where file are read from and write to
%
%   B1file-  the path to B1 file in epi space (this file was saved by
%            mrQ_smmothL_B1)
%
% OutPuts:
%
% AnalysisInfo -  a structure that keep records of the T1 fit files. the
%                 the ants registration files name and date will be recorded in it
%
%
%
% see also: mrQfit_T1M0_ver2 and mrQ_smmothL_B1.m



%% we first need make matlab work with the bachrc we have this is
% unforchanate complexisty


orig_path = getenv('LD_LIBRARY_PATH');

colon_idx = strfind(orig_path, ':');

matlab_path = [];
other_path = [];

for i=1:(length(colon_idx)-1)

    this_path = orig_path(colon_idx(i):colon_idx(i+1));
    % This is part of the matlab path:
    if strfind(this_path, 'matlab')
        matlab_path = [matlab_path, ':', this_path];

        % Should go before the matlab bit of the path (all the rest):
    else
        other_path = [other_path, ':', this_path];

    end
end

new_path = [other_path, ':', matlab_path];

setenv('LD_LIBRARY_PATH', new_path)

%%
 %% run ANTS non linear parameter to register SPGR to epi

%name the saved parameter file
AnalysisInfo.RB_WARP_EPI_SPGR=fullfile(outDir,'WARP_EPI_SPGR_RB');
%AnalysisInfo.RBMI_WARP_EPI_SPGR=fullfile(outDir,'WARP_EPI_SPGR_RBMI');


%if ~(exist([AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz'],'file') && exist([AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'],'file') )

    % ANTS: make the Ants comand
    cmANTS = ['xterm -e ANTS 3 -m CC[' t1fileHM ',' AnalysisInfo.T1_spgr_epi_RB ',1,2]  -o ' AnalysisInfo.RB_WARP_EPI_SPGR '.nii.gz --rigid-affine true'];
    %cmANTS=['xterm -e ANTS 3 -m MI[' t1fileHM ',' AnalysisInfo.T1_spgr_epi ',1,32]   -o ' AnalysisInfo.RBMI_WARP_EPI_SPGR '.nii.gz --rigid-affine true']
    % Run the command in unix and get back status and results:
    [status, result] = system(cmANTS,'-echo');
    if status ~=0
        cmANTS=['ANTS 3 -m CC[' t1fileHM ',' AnalysisInfo.T1_spgr_epi_RB ',1,2]  -o ' AnalysisInfo.RB_WARP_EPI_SPGR '.nii.gz --rigid-affine true'];
        [status, result] = system(cmANTS,'-echo')
    end
    


%%    name and saved sprg B1 map in spgr space
%the t1 file
AnalysisInfo.RB_B1_epi_spgr=fullfile(outDir,'Warp_B1_EPI2SPGR_RB.nii.gz');
%cmWarp=['xterm -e WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi ' '  AnalysisInfo.B1_epi_spgr ' -R '  B1file  ' ' AnalysisInfo.WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.WARP_EPI_SPGR 'Affine.txt'];

    % ANTS: make the Ants command
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' B1file ' '  AnalysisInfo.RB_B1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~=0
        cmWarp=['WarpImageMultiTransform  3 ' B1file ' '  AnalysisInfo.RB_B1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
        % Run the command in unix and get back status and results:
        [status, result] = system(cmWarp);
    end

%the B1 file
AnalysisInfo.RB_T1_epi_spgr=fullfile(outDir,'Warp_T1_EPI2SPGR_RB.nii.gz');

    % ANTS: make the Ants command
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi_RB ' '  AnalysisInfo.RB_T1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~= 0
        cmWarp=['WarpImageMultiTransform  3 ' AnalysisInfo.T1_spgr_epi_RB ' '  AnalysisInfo.RB_T1_epi_spgr ' -R '  t1fileHM  ' ' AnalysisInfo.RB_WARP_EPI_SPGR 'Warp.nii.gz ' AnalysisInfo.RB_WARP_EPI_SPGR 'Affine.txt'];
        [status, result] = system(cmWarp,'-echo');
    end


%%  % Reset the ld_path:
setenv('LD_LIBRARY_PATH', orig_path)

%% saved the update the analysis info file
infofile=fullfile(outDir,'AnalysisInfo.mat');
save(infofile,'AnalysisInfo');
