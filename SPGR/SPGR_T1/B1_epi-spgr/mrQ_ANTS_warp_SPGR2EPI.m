function [AnalysisInfo, Res]=mrQ_ANTS_warp_SPGR2EPI(SET1file,flipAngles,outDir,AlignFile,AligndSPGR,AnalysisInfo)
%
%mrQ_ANTS_warp_SPGR2EPI(AnalysisInfo,SET1file,flipAngles,outDir,AlignFile,AnalysisInfo)
% DOC!!!!
% NOTE: This function is currently misnamed. It is *not* a rigid body (RB)
% transformation.
%
%  This function evaluates the ANTs affine linear registration from SPGR to
%EPI. The ANTs will register the SPGR T1 to the SEIR EPI T1 using unix
%functions. The code assumes ANTs and ITK code are part of the computer
%unix sell path (./bashrc). The registration happens in two steps: First,
%define parameters; and second, apply it and warp the images. We warp both
%the T1 and the different flip angle degree SPGR raw images (all need to be
%aligned first).
%
% The process is documented in the AnalysisInfo structure.
%
% INPUTS:
%
%   SET1file     - The path from the SEIR T1 file to the SPGR will be 
%                  registered to this target.
%
%
%   flipAngles   - The raw data flip angles. (The raw data are saved as NIfTI 
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
%                  the ANTs registration files' name and date will be 
%                  recorded in it.
%
%   Res          - A structure of registered files, ready for B1 fit. This
%                  structure is also saved as the AlignFile.
%
% see also: mrQfit_T1M0_ver2
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016


 %% II. Run ANTs non-linear parameter to register SPGR to EPI


      
 %% IV. Name and save raw SPGR flipAngles images in EPI space  
           
 for d=1:length(flipAngles) %loop over flip Angles raw images
     
     % the name of the file we will align
     if (~exist('AligndSPGR','var')|| isempty(AligndSPGR)),
         rawIm=fullfile(outDir,['Align' num2str(flipAngles(d)) 'deg.nii.gz']);
     else
         rawIm=([AligndSPGR{d} '.nii.gz']) ;
         
     end
     
     %record it
     %  AnalysisInfo.Raw_spgr_epi{d}=fullfile(outDir,['WarpRB_' num2str(flipAngles(d)) 'deg_SPGRT2EPI.nii.gz']);
     AnalysisInfo.Raw_spgr_epi{d}=[AligndSPGR{d} 'WarpRB_SPGRT2EPI.nii.gz'];
     
     %ANTs: make the ANTs command
     cmWarp=['xterm -e WarpImageMultiTransform  3 ' rawIm  ' ' AnalysisInfo.Raw_spgr_epi{d} ' -R '  SET1file ' ' AnalysisInfo.WARP_SPGR_EPI 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI 'Affine.txt'];
     % Run the command in unix and get back status and results:
     [status, result] = system(cmWarp);
     if status ~= 0 ;
         cmWarp=['WarpImageMultiTransform  3 ' rawIm  ' ' AnalysisInfo.Raw_spgr_epi{d} ' -R '  SET1file ' ' AnalysisInfo.WARP_SPGR_EPI 'Warp.nii.gz ' AnalysisInfo.WARP_SPGR_EPI 'Affine.txt'];
         % Run the command in unix and get back status and results:
         [status, result] = system(cmWarp,'-echo')
     end
 end
        
%% V. Save and document 

%make a structure to work with for B1 fit
t1seir=readFileNifti(SET1file);
Res{1}.im=t1seir.data;
Res{1}.name='t1SEIRepi';
Res{1}.xform=t1seir.qto_xyz;
clear t1seir
t1spgr=readFileNifti(AnalysisInfo.T1_spgr_epi);
Res{2}.im=t1spgr.data;
Res{2}.name='t1SPGR_in_epi_space';

for i=1:length(flipAngles)
  im=readFileNifti(AnalysisInfo.Raw_spgr_epi{i}); 
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

