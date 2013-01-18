function location=mrQ_Callfs_autosegment(subjID, t1)
% Autosegment a t1 weighted anatomical scan using freesurfer
% autosegmentation, This function requires that you
% have freesurfer in your shell path.
%
% fs_autosegmentToITK(subjID, t1, [skipRecon], [resample_type])
%
% INPUTS:
%   subjID: directory name in which freesurfer stores outputs
%   t1:     file name (with complete path) of t1 used for segmentation.
%                TODO: Currently a single NIFTI is expected by this script.
%                Freesurfer is more flexible so the scipt could be improved
%                by allowing multiplte T1s or one or more directories of
%                DICOMs.
% 
% Example:
%   subjID  = 'jw';
%   t1      = fullfile('biac2/wandell2/data/anatomy/winawer/Anatomy20110308', 't1.nii.gz');
%   fs_autosegmentToITK(subjID, t1)
%
% 3/11/2009 Written by JW and HH.
% 6/30/2011 JW: added  'skipRecon' and 'resample_type' as optional input
%               arguments. Change the resample procedure to read the
%               deseried voxel size from the t1 file.
% 7/24/2012  AM:  edit for mrQ needs it only the call for freesurfer and it
%               in a difrent unix shell
%
%
% see also fs_ribbon2itk.m  fs_autosegmentToITK.m


%% Check Inputs & Paths

% subjID is the name of the directory that will be created by freesurfer to
% store segmentation and associated files.
if ~exist('subjID', 'var')
    warning('Subject ID is required input'); %#ok<WNTAG>
    eval('help fs_autosegmentToITK');
    return
end

% Get the path to the t1 file if it is not inputed. This can be any
% resolution. Freesurfer will resample to 1x1x1 for autosegmentation. We
% will extract a segmentation at the resolution of the original t1 at the
% end of this function.
if notDefined('t1') || ~exist(t1, 'file'),
    [fname pth] = uigetfile({'t1*.nii.gz', 'T1 files'; '*.nii', '.nii files'; '*.gz', '.gz files'}, 'Cannot locate T1 file. Please find it yourself.', pwd);
    t1 = fullfile(pth, fname);
end
if ~exist(t1, 'file'), error('Cannot locate t1 file'); end


% This is the directory where freesurfer puts subject data. If it is
% defined in the linux shell (e.g., bashrc) then matlab can find it. If it
% is not defined, look for the 'freesurfer_home/subjects', which is the
% default location in freesurfer.
subdir   = getenv('SUBJECTS_DIR');
if isempty(subdir),
    fshome = getenv('FREESURFER_HOME');
    subdir = fullfile(fshome, 'subjects');
end


location=fullfile(subdir,subjID) ;
if exist(location, 'dir')
    subjID=[subjID date]
location=fullfile(subdir,subjID) ;
end    
%% recon all (freesurfer will resample to 1 mm isotropic)
%
    fscmd = sprintf('xterm -e recon-all -i %s -subjid ''%s'' -all &', t1, subjID );
    
    [status result] =  system(fscmd);
     if status ~= 0
        fprintf('Something went wrong with the freesurfer call.');
        disp(result);
     end
%     
    
  %  eval(msg)
%%

