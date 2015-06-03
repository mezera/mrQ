function [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = ...
    mrQ_initSPGR_ver2(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ)
%
% [s,xform,mmPerVox] = ...
%   mrQ_initSPGR(spgrDir,[refImg],[mmPerVox],[interp],[skip],[clobber=false])
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
%       [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = mrQ_initSPGR(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ);
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
    spgrDir = uigetdir(pwd,'Select your SPGR directory');
end


% Rerference Image
if notDefined('refImg') || ~exist(refImg,'file')
    %   fprintf('\nNo reference image selected. ACPC alignemnt will be performed on raw data...\n')
elseif exist(refImg,'file')
    fprintf('Volumes will be aligned to: %s\n',refImg);
end

% Interpolation method
if notDefined('interp')
    interp = 1;
    mrQ.SPGR_init_interp=interp;
end
if interp ~= 1 && interp ~= 7
    error('Invalid interpolation method. Must be 1 (linear) or 7 (trilinear)');
end

% Clobber flag. Overwrite existing dat_aligned.mat if it already exists and
% redo the acpc.
if notDefined('clobber')
    clobber = false;
end

% Check if the acpc alignment should be done automatically or manually
if exist('mrQ','var') && ~isempty(mrQ) && isfield(mrQ,'autoacpc')
    autoAcpc = mrQ.autoacpc;
else
    autoAcpc = 0;
end


%
%%  Get paths to nifti SPGR data (get niiFiles)
% Get a structure that has the paths to each of the spgr nifti files.
if isfield(mrQ,'inputdata_spgr') %infut of list of nifti file and the relevant scan parameters
    [s niiFiles]=mrQ_input2Stuck(mrQ.inputdata_spgr,0);
    
    t1Inds(1:length(s))=1;
    t1Inds=logical(t1Inds);
    for ii = 1:numel(niiFiles)
        %  Load data from niftis - reshape and permute data (nifti) if needed
        
        [s1(ii)]= makeStructFromNifti(niiFiles{ii},-2,s(ii),mrQ.permution);
    end
    clear s
    s = s1;
    clear s1
    
else
    error('fdfdfdfd;')
end


%% More checks and housekeeping

% Get the flip angle for each of the series and store them in 'flip'. This
% will be passed into mrQ_multicoil_Weights *** Make sure this will
% work in the case that there are no nifti files (it should) ***
% We need this even if not processing the data again.




% Get the TR for the T1s [1 x numel(d)]
tr = [s(t1Inds).TR];


% Make sure that the TRs are all the same
if ~all(tr == tr(1))
    p = input(['TR''s do not match: ' num2str(tr) ' please select the TR, press 0 to use all of the scans  \n ']) ;
    if p==0
    else
        trIndsNot = [s(:).TR] ~= p;
        t1Inds    = t1Inds & ~trIndsNot;
    end
end


% Get the echo time for the t1s
te = [s(t1Inds).TE];


% Check that the TEs match
if ~all(te == te(1))
    p = input(['TE''s do not match: ' num2str(te) ' please select the TE, press 0 to use all of the scans  \n']) ;
    if p==0
    else
        teIndsNot = [s(:).TE] ~= p;
        t1Inds = t1Inds & ~teIndsNot;
    end
end


if ~isfield(mrQ,'cheack')
    mrQ.cheack=0;
end

if exist('skip','var')
    t1Inds(skip) = 0;
end

if mrQ.cheack==1
    for f = 1:numel(s(t1Inds)),
        showMontage(s(f).imData,[],[],[],[],10);
        
        an1 = input( ['Does the image  with of ' num2str(s(f).flipAngle) ' is good? Press 1 if yes 0 if no \n'])
        
        if an1==0
            t1Inds(f) = 0;
        end
        close figure 10
    end
    numData=length(find(t1Inds));
    
    if length(s)~=numData
        an = questdlg(['Do you like to continue the process? Are ' num2str(numData) ' scans enough data  ?'],' continue process ','YES','NO','YES');
        if strcmp('NO',an),
            error('the user stop the process');
        end
        mrQ.SPGR_Scan_Skiped_Num=find(t1Inds==0);
    end
    
end

%for f = 1:numel(s(t1Inds)), flipAngles(f) = s(f).flipAngle; end
niiFiles=niiFiles(t1Inds);
mrQ.SPGR_niiFile=niiFiles;
flipAngles=[s(t1Inds).flipAngle];
mrQ.SPGR_niiFile_FA=flipAngles;
% Set the resolution if it was not passed in
if ~exist('mmPerVox','var') || isempty(mmPerVox),
    mmPerVox = s(min(find(t1Inds))).mmPerVox(1:3);   %#ok<MXFND>
    mrQ.SPGR_init_mmPerVox=mmPerVox;
    
end
%end
%% ACPC Alignement

outDir = spgrDir;
% If the reference image was not passed in then we make the user take one
% of the images and choose the acpc landmarks and use the resulting image
% as a refernce for alignment.
if ~exist('refImg','var') || isempty(refImg)
    
    for i=1:numel(s),
        val(i) = 20-s(i).flipAngle; 
    end
    
    
    % Take a flip angle (closest to 20 deg) create a nifti image of that
    % volume and use for acpc marking
    [~, sec]     = sort(abs(val));
    fileRaw      = fullfile(outDir,'t1w_raw.nii.gz');
    t1w_acpcfile = fullfile(outDir,'t1w_acpc.nii.gz');
    dtiWriteNiftiWrapper(single(s(sec(1)).imData), s(sec(1)).imToScanXform, fileRaw);
    
    % Decide whether to do manual or automatic alignment
    if ~exist('autoAcpc','var') || isempty(autoAcpc) || autoAcpc == 0;
        % Do the acpc alignment - prompt the user to make sure it's good.
        an = 0;
        while an ~= 1 || isempty(an)
            mrAnatAverageAcpcNifti({fileRaw},t1w_acpcfile);
            an = questdlg('Does the alignment look good?','ACPC ALIGNMENT','YES','NO','YES');
            if strcmp('YES',an), an = 1; else an = 0; end
        end
    else
        % Automatically identify the ac and pc and mid sage by
        % computing a spatial normalization
        ni = readFileNifti(fileRaw);
        ni = niftiApplyCannonicalXform(ni);
        template =  fullfile(mrDiffusionDir, 'templates', 'MNI_T1.nii.gz');
        sn = mrAnatComputeSpmSpatialNorm(ni.data, ni.qto_xyz, template);
        c = mrAnatGetImageCoordsFromSn(sn, tal2mni([0,0,0; 0,-16,0; 0,-8,40])', true)';
        mrAnatAverageAcpcNifti(ni, t1w_acpcfile, c, [], [], [], false);
        
       
    end
    close all
    % The refImg is now the acpc aligned image.
    refImg = t1w_acpcfile;
    mrQ.SPGR_init_ref_acpc=refImg;
    
    % [s,xform] = relaxAlignAll(s(find(t1Inds)),[],mmPerVox,false,interp); *** WHAT'S THIS ***
end

%% ALIGNMENT: Do the alignment of the SPGRs and save out the aligned data

% Setup the output directory for the aligned data and make it if ~exist
outDir = fullfile(spgrDir,['Align_'  num2str(mmPerVox(1)) '_' num2str(mmPerVox(2)) '_'  num2str(mmPerVox(3))]);
if(~exist(outDir,'dir')), mkdir(outDir); end


% Align all the series to this subject's reference volume
%% NOTES
% originatly we used spm 8 ridge body registration code applayied
% by RFD for detail see relaxAlignAll. this was tested and work on
% GE data.
% this code is not working for our siemens data we there for use
% fsl implemntaion for details see mrQ_fslAlignCall. it still to test if all data need to be
% use this code.
%



% if  ~isfield(mrQ,'SPGR_raw_strac')
%spm
ref       = readFileNifti(refImg);

[s,xform] = relaxAlignAll(s(find(t1Inds)),ref,mmPerVox,true,interp); %#ok<FNDSB>
%else
%         Dpath=fullfile(mrQ.SPGR,'data');
% make nii from the images
%for ii=1:length(s)

%          niilist(ii).name=['T1w_FA' num2str(s(ii).flipAngle) '.nii.gz'];
%        savename=fullfile(Dpath,niilist(ii).name);
%  mrQ_makeNiftiFromStruct(s(ii),savename,s(ii).imToScanXform);
%end
% [s, xform ]=mrQ_fslAlignCall(Dpath,s,niilist,refImg);
% end

% Save out the aligned data
outFile = fullfile(outDir,'dat_aligned.mat');
save(outFile,'s', 'xform', 'mmPerVox');


mrQ.spgr_initDir=outDir;
mrQ.outDir=outDir;



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
