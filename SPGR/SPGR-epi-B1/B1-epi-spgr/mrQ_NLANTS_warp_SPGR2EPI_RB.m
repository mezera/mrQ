function [AnalysisInfo, Res]=mrQ_NLANTS_warp_SPGR2EPI_RB(SET1file,t1fileHM,flipAngles,outDir,AlignFile,AligndSPGR)
%
%mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)
%
% # this function evaluate the ANTS affine linear reistration from spgr to epi.
%the ANTS will register the spgr T1 to the SEIR epi T1 using unixs functions. the code esume ANTS and ITK code are part of the computer unix sell path (./bachrc).
%the registration happen in two steps first define parameters and secound
%apply it and wrap the images. we Warps both the T1 and the diferent flip angle degree SPGR raw  images(all need to be alignmened first).
% the procces is document in the  AnalysisInfo structure. 
%
% INPUTS:
%
%    AnalysisInfo -  a structure that keep records of the T1 fit files
%
%   SET1file  - the path fro  the SEIR T1 file the SPGR will be register to
%               this target
%
%   t1fileHM  - the path to SPGR T1 file (this is uncorected fast fit T1 map). this
%               file use to register to the SEIR target becouse they have similar
%                contrast (T1 maps).
%
%   flipAngles   %the raw data filp angles (the raw data are saved as nifft so
%                 scan parameter are not saved in there file but it is needed here)
%
%   outDir - the path to where file are read from and write to
%
%   AlignFile- the name of the saved output file 
%
%   AligndSPGR -the SPGR align Nifti files names
%
% OutPuts:
%
% AnalysisInfo -  a structure that keep records of the T1 fit files. the
%                 the ants registration files name and date will be recorded in it
%
% Res          -  a structure of register file ready for B1 fit. this
% structure is also saved as the AlignFile
%
%
% see also: mrQfit_T1M0_ver2



%% we first need make matlab work with the bachrc we have this is
% this was found to be unnecessary, but if ANTS doesn't seem to work -  
% uncomment this section to add it to the path, don't forget to uncomment 
% the path reset at the end of the scroipt:
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

 %% run ANTS non linear parameter to register SPGR to epi

%name the saved parameter file
AnalysisInfo.WARP_SPGR_EPI_RB= fullfile(outDir,'WARP_SPGR_EPI_RB');
 disp('ATNS  linear registration and Warp may take about 20 min ...')
 disp('...')
 
    % ANTS: make the Ants comand
    %cmANTS=['xterm -e ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -r Gauss[2,0] -o ' AnalysisInfo.WARP_SPGR_EPI '.nii.gz -i 30x99x11 -t SyN[0.5]'];
    cmANTS=['xterm -e ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -o ' AnalysisInfo.WARP_SPGR_EPI_RB '.nii.gz --rigid-affine true'];

     % Run the command in unix and get back status and results: 
    [status, result] = system(cmANTS);
    if status ~= 0 ; 
        cmANTS=['ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -o ' AnalysisInfo.WARP_SPGR_EPI_RB '.nii.gz --rigid-affine true'];
        [status, result] = system(cmANTS,'-echo')
    end


       
%%    name and saved sprg T1 map in epi space

% the T1 file name
AnalysisInfo.T1_spgr_epi_RB=fullfile(outDir,'Warp_T1_SPGRT2EPI_RB.nii.gz');

    % ANTS: make the Ants comand
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' t1fileHM ' ' AnalysisInfo.T1_spgr_epi_RB ' -R '  SET1file  ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];

    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~=0
        cmWarp=['WarpImageMultiTransform  3 ' t1fileHM ' ' AnalysisInfo.T1_spgr_epi_RB ' -R '  SET1file  ' ' AnalysisInfo.WARP_SPGR_EPI_RB 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI_RB 'Affine.txt'];
        [status, result] = system(cmWarp,'-echo')
    end
        
           %%    name and saved raw spgr flipAngles images in epi space  
           
           for d=1:length(flipAngles) %loof over  flip Angles raw images
               
               % the  name of the file we fill algin
               if (~exist('AligndSPGR','var')|| isempty(AligndSPGR)),
               rawIm=fullfile(outDir,['Align' num2str(flipAngles(d)) 'deg.nii.gz']);
               else
                                  rawIm=([AligndSPGR{d} '.nii.gz']) ;

               end
               
               %record it
             %  AnalysisInfo.Raw_spgr_epi_RB{d}=fullfile(outDir,['WarpRB_' num2str(flipAngles(d)) 'deg_SPGRT2EPI.nii.gz']);
               AnalysisInfo.Raw_spgr_epi_RB{d}=[AligndSPGR{d} 'WarpRB_SPGRT2EPI.nii.gz'];

               %ANTS: make the Ants comand
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

    
%make a structuure to work with for B1 fit
t1seir=readFileNifti(SET1file);
Res{1}.im=t1seir.data;
Res{1}.name='t1SEIRepi';
clear t1seir
t1spgr=readFileNifti(AnalysisInfo.T1_spgr_epi_RB);
Res{2}.im=t1spgr.data;
Res{2}.name='t1SPGR_in_epi_space';

for i=1:length(flipAngles)
  im=readFileNifti(AnalysisInfo.Raw_spgr_epi_RB{i}); 
Res{i+2}.im=im.data;
Res{i+2}.name=['align_rawFA' num2str(flipAngles(i))] ;

end;


%save the strucute
infofile=fullfile(outDir,'AnalysisInfo.mat');
save(infofile,'AnalysisInfo');
    save(AlignFile,'Res');
%     
%     % Reset the ld_path (unnecessary, see explanation at top of function)
% 
% setenv('LD_LIBRARY_PATH', orig_path)
% %saved the update the analysis info file


   %document  and done
AnalysisInfo.spgr2epi_Align_date=date;
disp('DONE with  linear registration SPGR to epi ')

