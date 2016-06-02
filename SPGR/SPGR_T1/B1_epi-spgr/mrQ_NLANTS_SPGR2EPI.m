function [WARP_image,WarpFiles]= mrQ_NLANTS_SPGR2EPI(LowResIm,HighResIm,maskfilepath,outDir,morefiles2_HR)%flipAngles,AlignFile,AligndSPGR)
%
%mrQ_NLANTS_SPGR2EPI(LowResIm,HighResIm,mask,flipAngles,outDir,AlignFile,AligndSPGR)
% !!!!!DOC
%
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
%   LowResIm     - The path to target for example the from the SEIR T1 file to the SPGR will be
%                  registered to this target.
%
%   HighResIm     - The path to high resulotion image. for example the SPGR T1 file (this is uncorrected fast-fit
%                  T1 map). This file is used to register to the SEIR target
%                  because they have similar contrast (T1 maps).
%   maskfilepath    a mask in low resultion space to register at.
%
%
%   outDir       - The path to which files are read from and written to.
%
%
%   morefiles2_HR   - The path to the hi-resultion files names to warp
%
%
% OUTPUTS:
%

%

%% I. We first need to make matlab work with the bashrc
% This was found to be unnecessary, but if ANTs doesn't seem to work,
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

%% II. Run ANTs non-linear parameter to register SPGR to EPI

%name the saved parameter file
WARP_image= fullfile(outDir,'MoveIm');
disp('ANTs linear registration and warp may take about 20 min ...')
disp('...')

% ANTS: make the ANTS command
%cmANTS=['xterm -e ANTS 3 -m CC[' SET1file ',' t1fileHM ',1,2] -r Gauss[2,0] -o ' AnalysisInfo.WARP_SPGR_EPI '.nii.gz -i 30x99x11 -t SyN[0.5]'];
if notDefined('maskfilepath')
    
    cmANTS=['xterm -e ANTS 3 -m CC[' LowResIm ',' HighResIm ',1,2] -o ' WARP_image '.nii.gz --rigid-affine true'];
    
else
    cmANTS=['xterm -e ANTS 3 -m CC[' LowResIm ',' HighResIm ',1,2] -x ' maskfilepath  ' -o ' WARP_image '.nii.gz --rigid-affine true'];
    % cmANTS=['xterm -e ANTS 3 -m CC[' maskfilepath ',' t1fileHM ',1,2] -o ' AnalysisInfo.WARP_SPGR_EPI_RB '.nii.gz --rigid-affine true'];
end
% Run the command in unix and get back status and results:
[status, result] = system(cmANTS);
if status ~= 0 ;
    if notDefined('maskfilepath')
        cmANTS=['ANTS 3 -m CC[' LowResIm ',' HighResIm ',1,2] -o ' WARP_image '.nii.gz --rigid-affine true'];
    else
        cmANTS=['ANTS 3 -m CC[' LowResIm ',' HighResIm ',1,2] -x ' maskfilepath  ' -o ' WARP_image '.nii.gz --rigid-affine true'];
        
    end
    [status, result] = system(cmANTS,'-echo')
end




%%    name and saved raw spgr flipAngles images in epi space
if ~notDefined('morefiles2_HR') 

for d=1:length(morefiles2_HR) %loof over  flip Angles raw images
    file=dir(morefiles2_HR{d});
    savefileN=fullfile(outDir,['Warp' file.name]);
    
    %ANTS: make the Ants comand
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' morefiles2_HR{d}  ' ' savefileN ' -R '  LowResIm ' ' WARP_image 'Warp.nii.gz ' WARP_image 'Affine.txt'];
    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~= 0 ;
        cmWarp=['WarpImageMultiTransform  3 ' morefiles2_HR{d}  ' ' savefileN ' -R '  LowResIm ' ' WARP_image 'Warp.nii.gz ' WARP_image 'Affine.txt'];
        % Run the command in unix and get back status and results:
        [status, result] = system(cmWarp,'-echo')
    end
    
    WarpFiles{d}=savefileN;
    
end

end



