function [s,xform,mmPerVox,niiFiles,mrQ] = ...
    mrQ_initMT(spgrDir,refImg,mmPerVox,interp,skip,mtOffset,clobber,mrQ)
%
% [s,xform,mmPerVox] = ...
%   mrQ_initSPGR_ver2(spgrDir,[refImg],[mmPerVox],[interp],[skip],[clobber=false])
%
% Load and align all the SPGR dicoms in the spgrDir. (It won't work
% for multi-coil data, so take care of that in mrQ_arrangeData.m).
%
% This function will align dicoms with the same TR and TE (for T1-M0 linear
% fit). So if there is a different TR or TE you will have to select the the
% right TR and TE. Note that you can also work with multiple TR's but then
% the relevant line needs to be commented. *** MORE INFO ON THIS NEEDED ***
%.
%
%
% INPUTS
%       spgrDir:  Where the SPGR dicom directories are
%
%       refImg:   Different ref images can be used as an input (refImg is a
%                 path to a nifti image). If there is no refImg (refImg is
%                 empty) then the SPGR with a similar contrast to the T1
%                 weighted image will be selected and the user will be
%                 asked to mark the ac/pc using mrAnatAverageAcpcNifti
%
%       mmPerVox: The resolution at which you want to resample the data.
%                 This is a 3X1 (1x3?) vector. If empty, the dicom
%                 resolution will be used-this does not have to be the
%                 native scan size as the magnet output % is zeroed. The
%                 saved directory will have the resolution in its name.
%
%       interp:   Interpolation method. [Default = 1]
%                 1 = trilinear,
%                 7 = b-spline (resampling algorithm)
%
%       skip:     you can skip any of the scans in spgrDir if you want by
%                 passing in a 1xn vector of scans to skip.
%
%
%      clobber:   Overwrite existing data and reprocess. [default = false]
%
%        mrQ:     information structure

%
%
% OUTPUTS
%       S:        All the aligned images are saved within the structure S.
%
%       xform:    The matrix that transforms from raw to acpc.
%
%       mmPerVox: The resolution that the data was sampled to.
%
%                 * These three outputs are saved in the output directory
%                   in the file: fullfile(outDir,'dat_aligned.mat');
%
%       niiFiles: A list of full paths to the nifti files used. If nifti
%                 files were not used this will be = [];
%
%     flipAngles: A 1xn vector of flip angles.
%
%     mrQ:          information structure, updated

%
%
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
%
% EXAMPLE USAGE
%       spgrDir  = '/viridian/scr1/data/qmr/20111109_1426/SPGR_1';
%       refImg   = [];
%       mmPerVox = [];
%       interp   = [];
%       skip     = [];
%       clobber  = [];
%       [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = mrQ_initSPGR_ver2(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ);
%
%
% (C) Stanford University, VISTA Lab
%

%#ok<*NODEF>
%#ok<*AGROW>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To Do:
%if there are no niffti we need to make them from the dicom and not work with the dicoms after that
%    that wiil be more compatible with the rest of the function!!! in the
%    CNi we always have dicom so it not a big problem for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check INPUTS

% Data directory
if notDefined('spgrDir') || ~exist(spgrDir,'dir')
    spgrDir = uigetdir(pwd,'Select your MT directory');
end


% Rerference Image
if notDefined('refImg') || ~exist(refImg,'file')
    fprintf('\nNo reference image selected. ACPC alignemnt will be performed on raw data...\n')
elseif exist(refImg,'file')
    fprintf('Volumes will be aligned to: %s\n',refImg);
end

% Interpolation method
if notDefined('interp')
    interp = 1;
     mrQ.MT.MT_init_interp=interp;
end
if interp ~= 1 && interp ~= 7
    error('Invalid interpolation method. Must be 1 (linear) or 7 (trilinear)');
end




% Clobber flag. Overwrite existing dat_aligned.mat if it already exists and
% redo the acpc.
if notDefined('clobber')
    clobber = false;
end



 process = true;
%%  Get paths to nifti SPGR data (get niiFiles)
% Get a structure that has the paths to each of the spgr nifti files.

dicomDir = fullfile(spgrDir,'data');
niiFiles = [];


% Get the path to each of the dicom directories (d)
d = genpath(dicomDir);
if(isempty(d)), error(['Dicom dir "' d '" not found or empty.']); end
if(isunix), d = explode(':',genpath(dicomDir));
else d = explode(';',genpath(dicomDir)); end


% The first and last entries in d are the root and an empty matrix
% (resp), so we get rid of those so that d is now simply the
% directories contating the dicom/niftis.
d = d(2:end-1);


% Check to see is there is a nifti file in each of the SPGR
% directories. If there is then status = 0 - a bit hacky.
cmd = ['ls -R ' dicomDir '/*/*nii.gz'];
[status,~] = system(cmd);


% Nifti files are in the dicom directories - load from there.
if status == 0
    % Loop over 'd' to get the nifi paths DO THIS EARLIER - OUTSIDE OF THIS
    % STATEMENT
    for ii = 1:numel(d)
        ni = dir(fullfile(d{ii},'*nii.gz'));
        niiFiles{ii} = fullfile(d{ii},ni.name);
        % 'niiFiles' now contains all the paths to the nifti files and can
        % be passed into mrQ_multicoil_Weights
    end
    
else
    % This needs to look in the raw directory for the nifti files - try to
    % guess, prompt if not found.
    rawDir = fullfile(mrvDirup(spgrDir),'raw');
    if ~exist(rawDir,'dir')
        rawDir = uigetdir(spgrDir,'Set raw directory');
    end
    
    for ii = 1:numel(d)
        [~, f ~] = fileparts(d{ii});
        ni = dir(fullfile(rawDir,[f '.nii.gz']));
        niiFiles{ii} = fullfile(rawDir,ni.name);
        
    end
    
end

% a.m  is this a bug i don't think it works is it  needed?

% Make sure that niiFiles matches up with s(n).flipAngle
% Check that the numel are the same.
if ~process
    for nn=1:numel(s)
        fName = dir(fullfile(rawDir,['*' num2str(s(nn).flipAngle) 'deg*.nii.gz']));
        tmp{nn} = fullfile(rawDir,fName.name);
    end
    if ~isequal(niiFiles,tmp)
        disp('Path string order in niiFiles does not match that in s. Fixing.');
        niiFiles = tmp;
    end
end


%% Load data from niftis or dicoms - reshape and permute data (nifti)

if process
    
    % Make a dummy structure so we will have the dicom info
    s = dicomLoadAllSeries(dicomDir);
    
    % Loop over niiFiles to get the data from the nifti and combine
    % with the dicom info - reshape and permute.
    for ii = 1:numel(niiFiles)
        [s1(ii) mrQ.MT.coilNum(ii)]= makeStructFromNifti(niiFiles{ii},[],s(ii));
    end
    clear s
    s = s1;
    clear s1
end




mrQ.MT.MT_niiFile=niiFiles;


if process
    
    % Check which volumes match SPGR sequence name: '3DGRASS' or 'EFGRE3D'
%     spgrInds(1:length(s)) = 0;
%     sequenceNames = {s(:).sequenceName};
%     spgrInds(strcmp('3DGRASS',sequenceNames)) = 1;
%     spgrInds(strcmp('EFGRE3D',sequenceNames)) = 1;
%     
%     
%     % Find MT & T1 indices
%     mtInds = [s(:).mtOffset] > 0;
%     t1Inds = spgrInds & ~mtInds;
%     
    
    % Get the TR for the T1s [1 x numel(d)]
    tr = [s(:).TR];
    
    
    % Make sure that the TRs are all the same
%     if ~all(tr == tr(1))
%         p = input(['TR''s do not match: ' num2str(tr) ' please select the TR ']) ;
%         trIndsNot = [s(:).TR] ~= p;
%         t1Inds    = t1Inds & ~trIndsNot;
%     end
    
    
    % Get the echo time for the t1s
    te = [s(:).TE];
    
    
    % Check that the TEs match
%     if ~all(te == te(1))
%         p = input(['TE''s do not match: ' num2str(te) ' please select the TE ']) ;
%         teIndsNot = [s(:).TE] ~= p;
%         t1Inds = t1Inds & ~teIndsNot;
%     end
    
    
    % Skip scans
    if exist('skip','var')
        t1Inds(skip) = 0;
    end

    
    % Set the resolution if it was not passed in
    if ~exist('mmPerVox','var') || isempty(mmPerVox),
        mmPerVox = s(1).mmPerVox(1:3);   % #ok<MXFND>
            mrQ.MT.MT_init_mmPerVox=mmPerVox;

    end
end
%% ACPC Alignement
if process
    outDir = spgrDir;
    % If the reference image was not passed in then we make the user take one
    % of the images and choose the acpc landmarks and use the resulting image
    % as a refernce for alignment.
    if ~exist('refImg','var') || isempty(refImg)
        
        f
        fileRaw      = fullfile(outDir,'t1w_raw.nii.gz');
        t1w_acpcfile = fullfile(outDir,'t1w_acpc.nii.gz');
        dtiWriteNiftiWrapper(single(s(1).imData), s(1).imToScanXform, fileRaw);
        
        
        % Do the acpc alignment - prompt the user to make sure it's good.
        an = 0;
        while an ~= 1 || isempty(an)
            mrAnatAverageAcpcNifti({fileRaw},t1w_acpcfile);
            an = questdlg('Does the alignment look good?','ACPC ALIGNMENT','YES','NO','YES');
            if strcmp('YES',an), an = 1; else an = 0; end
        end
        
        close all
        % The refImg is now the acpc aligned image.
        refImg = t1w_acpcfile;
        mrQ.MT.MT_init_ref_acpc=refIm;

        % [s,xform] = relaxAlignAll(s(find(t1Inds)),[],mmPerVox,false,interp); *** WHAT'S THIS ***
    end
end


for i=1:length(mtOffset)
    s(i).mtOffset=mtOffset(i);
end

%% ALIGNMENT: Do the alignment of the SPGRs and save out the aligned data

if process
    % Setup the output directory for the aligned data and make it if ~exist
    outDir = fullfile(spgrDir,['Align_'  num2str(mmPerVox(1)) '_' num2str(mmPerVox(2)) '_'  num2str(mmPerVox(3))]);
    if(~exist(outDir,'dir')), mkdir(outDir); end
    
    
    % Align all the series to this subject's reference volume
    ref       = readFileNifti(refImg);
    [s,xform] = relaxAlignAll(s(:),ref,mmPerVox,true,interp); %#ok<FNDSB>
    
    
    % Save out the aligned data
    outFile = fullfile(outDir,'dat_aligned.mat');
    save(outFile,'s', 'xform', 'mmPerVox');
    
    
    mrQ.MT.MT_initDir=outDir;
       

end


%% Sanity check: make sure that the order in niiFiles matches s.flipAngles

% Make sure that niiFiles matches up with s(n).flipAngle
% Check that the numel are the same.
% fprintf('Verifying flipAngle ordering ...');
%    
% % i deside not to use this cheack for now
% load (outFile)
% 
% for nn=1:numel(s)
%     fName = dir(fullfile(mrvDirup(niiFiles{1}),['*' num2str(s(nn).flipAngle) 'deg*.nii.gz']));
%     temp{nn} = s(nn).flipAngle;
% end
% if ~isequal(niiFiles,temp)
%     disp('Path string order in niiFiles does not match that in s. Fixing.');
%     niiFiles = temp;
% end
% fprintf('Done.\n');






return
