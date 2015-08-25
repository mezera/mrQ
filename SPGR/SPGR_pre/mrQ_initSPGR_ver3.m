function [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = ...
    mrQ_initSPGR_ver3(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ)
%
% [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = ...
%     mrQ_initSPGR_ver2(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ)
% 
% Loads and aligns all the SPGR NIfTIs in the spgrDir. (It won't work
% for multi-coil data, so take care of that in mrQ_arrangeData.m).
%
% This function will align NIfTIs with the same TR and TE (for T1-M0 linear
% fit). So if there is a different TR or TE, you will have to select the
% correct TR and TE. Note that you can also work with multiple TR's, but
% then the relevant line needs to be commented. A nonlinear function can
% account for different TR's and/or TE's, but this function does not do
% that, though it may be a future direction for additional coding.
%
%
% INPUTS
%       spgrDir:  The SPGR NIfTI output directory
%
%       refImg:   Different ref images can be used as an input (refImg is a
%                 path to a NIfTI image). If there is no refImg (i.e. refImg 
%                 is empty), then the SPGR with a similar contrast to the
%                 T1-weighted image will be selected and the user will be
%                 asked to mark the ac-pc using mrAnatAverageAcpcNifti.m.
%
%       mmPerVox: The resolution at which you want to resample the data.
%                 This is a 3X1 vector. If empty, the NIfTI
%                 resolution will be used -- this does not have to be the
%                 native scan size, as the magnet output is zeroed. The
%                 saved directory will have the resolution in its name.
%
%       interp:   Interpolation method.
%                 1 = trilinear (default)
%                 7 = b-spline (resampling algorithm)
%
%       skip:     You can skip any of the scans in spgrDir by entering
%                 a 1xN vector of scans to skip.
%
%       clobber:  Overwrite existing data and reprocess. [Default = false]
%
%       mrQ:      Information structure.
%
%
%
% OUTPUTS
%     s:          All the aligned images are saved within the structure s.
%
%     xform:      The matrix that transforms from raw to ac-pc.
%
%     mmPerVox:   The resolution to which the data was sampled.
%
%                 * These three outputs are saved in the output directory
%                   in the file: fullfile(outDir,'dat_aligned.mat');
%
%     niiFiles:   A list of full paths to the NIfTI files used. If NIfTI
%                  files were not used, this will be = [];
%
%     flipAngles: A 1xN vector of flip angles.
%
%     mrQ:        Information structure, updated
%
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
%       [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] = ...
%                          mrQ_initSPGR_ver2(spgrDir,refImg,mmPerVox,interp,skip,clobber,mrQ);
%
%
% (C) Stanford University, VISTA Lab
%

%#ok<*NODEF>
%#ok<*AGROW>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To Do:
%    If there are no NIfTIs, we need to make them from the DICOMs. We will use
%    NIfTIs (and not DICOMs), as they are more compatible with the rest of
%    the function! In the CNI we always have DICOM so it is not a big
%    problem for now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I. Check INPUTS

% Data directory
if notDefined('spgrDir') || ~exist(spgrDir,'dir')
    spgrDir = uigetdir(pwd,'Select your SPGR directory');
end


% Reference Image
if notDefined('refImg') || ~exist(refImg,'file')
    %   fprintf('\n  No reference image selected. ACPC alignemnt will be performed on raw data...\n')
elseif exist(refImg,'file')
    fprintf('Volumes will be aligned to: %s\n',refImg);
end

% Interpolation method
if notDefined('interp')
    interp = 1;
    mrQ.SPGR_init_interp=interp;
end

if interp ~= 1 && interp ~= 7
    error('Invalid interpolation method. Must be 1 (trilinear) or 7 (b-spline)');
end

% Clobber flag. Overwrite existing dat_aligned.mat if it already exists and
% redo the ac-pc.
if notDefined('clobber')
    clobber = false;
end

% Check if the ac-pc alignment should be done automatically or manually
if exist('mrQ','var') && ~isempty(mrQ) && isfield(mrQ,'autoacpc')
    autoAcpc = mrQ.autoacpc;
else
    autoAcpc = 0;
end


%%  II. Get paths to NIfTI SPGR data (get niiFiles)
% Get a structure that has the paths to each of the SPGR NIfTI files.

if isfield(mrQ,'inputdata_spgr') %input of list of nifti files and the relevant scan parameters
    [s, niiFiles]=mrQ_input2struct(mrQ.inputdata_spgr,0);
    for ii=1:4; s(ii).imData=s(ii).imData(:,:,:,1);end;
    t1Inds(1:length(s))=1;
    t1Inds=logical(t1Inds);
    for ii = 1:numel(niiFiles)
        %  Load data from niftis - reshape and permute data (nifti) if needed
        
        [s1(ii)]= makeStructFromNifti(niiFiles{ii},-2,s(ii),mrQ.permutation);
    end
    clear s
    s = s1;
    clear s1
    
else
    error('mrQ.inputdata_spgr does not exist! Please define the list of NIfTI files and the relevant scan parameters. ')
end


%% III. More checks and housekeeping

% Get the flip angle for each of the series and store them in 'flip'. This
% will be passed into mrQ_multicoil_Weights *** Make sure this will
% work in the case that there are no nifti files (it should) ***
% We need this even if not processing the data again.


% Get the TR for the T1s [1 x numel(d)]
tr = [s(t1Inds).TR];

% Make sure that the TRs are all the same
if ~all(tr == tr(1))
    p = input(['TR''s do not match: ' num2str(tr) ' please select the TR, press 0 to use all of the scans  \n ']) ;
    if p~=0
        trIndsNot = [s(:).TR] ~= p;
        t1Inds    = t1Inds & ~trIndsNot;
    end
end


% Get the echo time for the t1s
te = [s(t1Inds).TE];


% Check that the TEs match
if ~all(te == te(1))
    p = input(['TE''s do not match: ' num2str(te) ' Please select the TE, press 0 to use all of the scans  \n']) ;
    if p~=0
        teIndsNot = [s(:).TE] ~= p;
        t1Inds = t1Inds & ~teIndsNot;
    end
end


if ~isfield(mrQ,'check')
    mrQ.check=0;
end

if exist('skip','var')
    t1Inds(skip) = 0;
end

if mrQ.check==1
    for f = 1:numel(s(t1Inds)),
        showMontage(s(f).imData,[],[],[],[],10);
        
        an1 = input( ['Is the image  with the flip angle of ' num2str(s(f).flipAngle) ' good? Press 1 if yes, 0 if no \n'])
        
        if an1==0
            t1Inds(f) = 0;
        end
        close figure 10
    end
    
    numData=length(find(t1Inds));
    
    if length(s)~=numData
        an = questdlg(['Would you like to continue the process? Are ' num2str(numData) ' scans enough data  ?'],' continue process ','YES','NO','YES');
        if strcmp('NO',an),
            error('The user stopped the process.');
        end
        mrQ.SPGR_Scan_Skiped_Num=find(t1Inds==0);
    end
    
end

%for f = 1:numel(s(t1Inds)), flipAngles(f) = s(f).flipAngle; end
niiFiles=niiFiles(t1Inds);
mrQ.SPGR_niiFile=niiFiles;
flipAngles=[s(t1Inds).flipAngle];
mrQ.SPGR_niiFile_FA=flipAngles;

te = [s(t1Inds).TE];
tr = [s(t1Inds).TR];
mrQ.SPGR_niiFile_TR=tr;
mrQ.SPGR_niiFile_TE=te;


% Set the resolution if it was not entered in
if ~exist('mmPerVox','var') || isempty(mmPerVox),
    mmPerVox = s(min(find(t1Inds))).mmPerVox(1:3);   %#ok<MXFND>
    mrQ.SPGR_init_mmPerVox=mmPerVox; 
end
%end

%% IV. AC-PC Alignment

outDir = spgrDir;
% If the reference image was not entered in, then we make the user take one
% of the images and choose the ac-pc landmarks. We then use the resulting image
% as a reference for alignment.

if ~exist('refImg','var') || isempty(refImg)
    
    for i=1:numel(s),
        val(i) = 20-s(i).flipAngle; 
    end
    
    
    % Take a flip angle (closest to 20 deg) create a NIfTI image of that
    % volume and use for ac-pc marking
    [~, sec]     = sort(abs(val));
    fileRaw      = fullfile(outDir,'t1w_raw.nii.gz');
    t1w_acpcfile = fullfile(outDir,'t1w_acpc.nii.gz');
    dtiWriteNiftiWrapper(single(s(sec(1)).imData), s(sec(1)).imToScanXform, fileRaw);
    
    % Decide whether to do manual or automatic alignment
    if ~exist('autoAcpc','var') || isempty(autoAcpc) || autoAcpc == 0;
        % Do the ac-pc alignment - prompt the user to make sure it's good.
        an = 0;
        while an ~= 1 || isempty(an)
            mrAnatAverageAcpcNifti({fileRaw},t1w_acpcfile);
            an = questdlg('Does the alignment look good?','ACPC ALIGNMENT','YES','NO','YES');
            if strcmp('YES',an), an = 1; else an = 0; end
        end
    else
        % Automatically identify the AC and PC midsagittally by
        % computing a spatial normalization
        ni = readFileNifti(fileRaw);
        ni = niftiApplyCannonicalXform(ni);
        template =  fullfile(mrDiffusionDir, 'templates', 'MNI_T1.nii.gz');
        sn = mrAnatComputeSpmSpatialNorm(ni.data, ni.qto_xyz, template);
        c = mrAnatGetImageCoordsFromSn(sn, tal2mni([0,0,0; 0,-16,0; 0,-8,40])', true)';
        mrAnatAverageAcpcNifti(ni, t1w_acpcfile, c, [], [], [], false);
        
       
    end
    
    close all
    
    % The refImg is now the ac-pc aligned image.
    refImg = t1w_acpcfile;
    mrQ.SPGR_init_ref_acpc=refImg;
    
    % [s,xform] = relaxAlignAll(s(find(t1Inds)),[],mmPerVox,false,interp); *** WHAT'S THIS ***
end

%% V. ALIGNMENT
% Do the alignment of the SPGRs and save out the aligned data

% Setup the output directory for the aligned data and make it if ~exist
outDir = fullfile(spgrDir,['Align_'  num2str(mmPerVox(1)) '_' num2str(mmPerVox(2)) '_'  num2str(mmPerVox(3))]);
if(~exist(outDir,'dir')), mkdir(outDir); end


% Align all the series to this subject's reference volume



%% NOTES
% Originally we used spm 8 rigid body registration code applied by RFD; for
% details, see relaxAlignAll. This was tested and it works on GE data but
% not for Siemens data. Therefore we use fsl implementation; for details, see
% mrQ_fslAlignCall. It still needs to be tested for other kinds of data.
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
%          savename=fullfile(Dpath,niilist(ii).name);
%          mrQ_makeNiftiFromStruct(s(ii),savename,s(ii).imToScanXform);
%end
% [s, xform ]=mrQ_fslAlignCall(Dpath,s,niilist,refImg);
% end

% Save out the aligned data
outFile = fullfile(outDir,'dat_aligned.mat');
save(outFile,'s', 'xform', 'mmPerVox');

 mrQ.spgr_initDir=outDir;

 mrQ.xform=xform;

% save the aligned images in the mrQ struct
for j=1:length(s)
    ref   = fullfile(outDir,['Align' num2str(flipAngles(j)) 'deg']);
    dtiWriteNiftiWrapper(single(s(j).imData), xform,ref);
    mrQ.AligndSPGR{j}=ref;
    
end


%% Sanity check: make sure that the order in niiFiles matches s.flipAngles

% Make sure that niiFiles match up with s(n).flipAngle
% Check that the numel are the same.
% fprintf('Verifying flipAngle ordering ...');
%
% % I decided not to use this check for now
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
