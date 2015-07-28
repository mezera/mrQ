function [AnalysisInfo, Res]=mrQ_NLANTS_warp_SPGR2EPI_RB(SET1file,t1fileHM,flipAngles,outDir,AlignFile,AligndSPGR)
%
%mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)
%
%  This function evaluates the ANTS affine linear registration from SPGR to
%EPI. The ANTS will register the SPGR T1 to the SEIR EPI T1 using unix
%functions. The code assumes ANTS and ITK code are part of the computer unix
%sell path (./bachrc). The registration happens in two steps: First, define
%parameters; and second, apply it and warp the images. We warp both the T1
%and the different flip angle degree SPGR raw images (all need to be
%aligned first).
%
% The process is documented in the AnalysisInfo structure.
%
% INPUTS:
%
%   SET1file     - The path from the SEIR T1 file to the SPGR will be 
%                  registered to this target.
%
%   t1fileHM     - The path to the SPGR T1 file (this is uncorected fast-fit
%                  T1 map). This file is used to register to the SEIR target 
%                  because they have similar contrast (T1 maps).
%
%   flipAngles   - The raw data flip angles (the raw data are saved as NIfTI 
%                  so the scan parameters are not saved in that file, 
%                  but it is needed here).
%
%   outDir       - The path to which files are read from and written to.
%
%   AlignFile    - The name of the saved output file. 
%
%   AligndSPGR   - The SPGR-aligned NIfTI files names.
%
%
% OUTPUTS:
%
%   AnalysisInfo - A structure that keeps records of the T1 fit files. The
%                  the ANTS registration files' name and date will be 
%                  recorded in it.
%
%   Res          - A structure of registered files, ready for B1 fit. This
%                  structure is also saved as the AlignFile.
%
% see also: mrQfit_T1M0_ver2
%

%% We first need to make matlab work with the bachrc
% This was found to be unnecessary, but if ANTS doesn't seem to work,  
% un-comment this section to add it to the path. Don't forget to un-comment 
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
%    this_path = orig_path(colon_idx(i):colon_idx(i+1));
%    % This is part of the matlab path:
%    if strfind(this_path, 'matlab')
%        matlab_path = [matlab_path, ':', this_path];
% 
%    % Should go before the matlab bit of the path (all the rest):
%    else
%        other_path = [other_path, ':', this_path];
% 
%    end
% end
% 
% new_path = [other_path, ':', matlab_path];
% 
% setenv('LD_LIBRARY_PATH', new_path)

 %% run ANTS non-linear parameter to register SPGR to EPI

%name the saved parameter file
AnalysisInfo.WARP_SPGR_EPI_RB= fullfile(outDir,'WARP_SPGR_EPI_RB');
 disp('ANTS linear registration and warp may take about 20 min ...')
 disp('...')
 
    % ANTS: make the ANTS command
    %cmANTS=['xterm -e ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -r Gauss[2,0] -o ' AnalysisInfo.WARP_SPGR_EPI '.nii.gz -i 30x99x11 -t SyN[0.5]'];
    cmANTS=['xterm -e ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -o ' AnalysisInfo.WARP_SPGR_EPI_RB '.nii.gz --rigid-affine true'];

     % Run the command in unix and get back status and results: 
    [status, result] = system(cmANTS);
    if status ~= 0 ; 
        cmANTS=['ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -o ' AnalysisInfo.WARP_SPGR_EPI_RB '.nii.gz --rigid-affine true'];
        [status, result] = system(cmANTS,'-echo')
    end


%%    name and save SPGR T1 map in EPI space

% the T1 file name
AnalysisInfo.T1_spgr_epi_RB=fullfile(outDir,'Warp_T1_SPGRT2EPI_RB.nii.gz');

    % ANTS: make the ANTS command
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' t1fileHM ' ' AnalysisInfo.T1_spgr_epi_RB ' -R '  SET1file  ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];

    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~=0
        cmWarp=['WarpImageMultiTransform  3 ' t1fileHM ' ' AnalysisInfo.T1_spgr_epi_RB ' -R '  SET1file  ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];
        [status, result] = system(cmWarp,'-echo')
    end
        
 %%    name and save raw SPGR flipAngles images in EPI space  
           
           for d=1:length(flipAngles) %loop over flip Angles raw images
               
               % the name of the file we fill align
               if (~exist('AligndSPGR','var')|| isempty(AligndSPGR)),
               rawIm=fullfile(outDir,['Align' num2str(flipAngles(d)) 'deg.nii.gz']);
               else
                                  rawIm=([AligndSPGR{d} '.nii.gz']) ;

               end
               
               %record it
             %  AnalysisInfo.Raw_spgr_epi_RB{d}=fullfile(outDir,['WarpRB_' num2str(flipAngles(d)) 'deg_SPGRT2EPI.nii.gz']);
               AnalysisInfo.Raw_spgr_epi_RB{d}=[AligndSPGR{d} 'WarpRB_SPGRT2EPI.nii.gz'];

               %ANTS: make the ANTS command
               cmWarp=['xterm -e WarpImageMultiTransform  3 ' rawIm  ' ' AnalysisInfo.Raw_spgr_epi_RB{d} ' -R '  SET1file ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];
               % Run the command in unix and get back status and results:
               [status, result] = system(cmWarp);
               if status ~= 0 ;
                   cmWarp=['WarpImageMultiTransform  3 ' rawIm  ' ' AnalysisInfo.Raw_spgr_epi_RB{d} ' -R '  SET1file ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];  
                   % Run the command in unix and get back status and results:
                   [status, result] = system(cmWarp,'-echo')
               end
           end
        
%% save and document 

    
%make a structure to work with for B1 fit
t1seir=readFileNifti(SET1file);
Res{1}.im=t1seir.data;
Res{1}.name='t1SEIRepi';
Res{1}.xform=t1seir.qto_xyz;
clear t1seir
t1spgr=readFileNifti(AnalysisInfo.T1_spgr_epi_RB);
Res{2}.im=t1spgr.data;
Res{2}.name='t1SPGR_in_epi_space';

for i=1:length(flipAngles)
  im=readFileNifti(AnalysisInfo.Raw_spgr_epi_RB{i}); 
Res{i+2}.im=im.data;
Res{i+2}.name=['align_rawFA' num2str(flipAngles(i))] ;

end;


%save the structure
    save(AlignFile,'Res');
%     
%     % Reset the ld_path (unnecessary, see explanation at top of function)
% 
% setenv('LD_LIBRARY_PATH', orig_path)
% %saved the update the analysis info file


   %document and done
AnalysisInfo.spgr2epi_Align_date=date;
disp('DONE with linear registration SPGR to EPI ')

