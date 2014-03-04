function [AnalysisInfo]=mrQ_CSF(outDir,freesurfer,T1file,AnalysisInfo)
%
% Create a CSF ROI from the freesurfer segmentaiton
% 
% mrQ_CSF(outDir,freesurfer,T1file)
%
% This function will generate a smooth white matter mask for for the coil
% gain estimatin with polynomial with a degree (defult 3) and save a clean
% PD and WF maps.
% 
% Provide a freesurfer raw segmentation image in nii.gz format
% (aparc+aseg.nii.gz ) and the directory of the data analysis and a M0file
% and T1 images.
%
% INPUTS:
%       outDir      - the directory of the data analysis
%       freesurfer  - Freesufarer raw segmentation image in nii.gz format
%                     (aparc+aseg.nii.gz )
%       T1File      - E.g., 'T1_lsq_GLr.nii.gz'
% 
% OUTPUTS:
%       Saves out 'FS_tissue.nii.gz' and 'csf_FS_T1.nii.gz' to the outDir.
% 
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%       
% 
% (C) Stanford University, VISTA Lab 
% 

%#ok<*FNDSB>


%% CHECK INPUTS

if notDefined('outDir')
    outDir = uigetdir(pwd,'Choose your analysis directory');
end

if exist('T1file','var') && ~isempty(T1file)
    disp(['Loading T1 data from ' T1file '...']);
    T1 = readFileNifti(T1file);
    T1 = double(T1.data);
else
    T1file=mrQ_getT1file(mrQ);
   % T1file = fullfile(outDir,'T1_map_lsq.nii.gz');
   % disp(['trying  to load T1 from ' T1file '...']);
    if(exist(T1file,'file'))
        disp(['Loading T1 data from ' T1file '...']);
        T1 = readFileNifti(T1file);
        T1 = double(T1.data);
    
   %else
    %    T1file = mrvSelectFile('r','*.nii.gz','Select T1 File');
        % disp(['error , can not find the file: '  T1file]);
        %error
    end
end

if notDefined('freesurfer') || ~exist(freesurfer,'file')
    freesurfer = mrvSelectFile('r','Select Freesurfer segmentation');
end

BMfile = fullfile(outDir,'brainMask.nii.gz');
if ~exist(BMfile,'file')
    BMfile = mrvSelectFile('r','Select the Brain Mask');
end

%% Load Data

% Load the Analysis info file and append it with the info provided
infofile = fullfile(outDir,'AnalysisInfo.mat');
load(infofile);

% AnalysisInfo.GainPolydegrees = degrees;
% AnalysisInfo.M0forWF = M0file;

AnalysisInfo.T1wSeg  = freesurfer;
AnalysisInfo.T1forWF        = T1file;
AnalysisInfo.WFdate         = date;

save(infofile,'AnalysisInfo');

% Load the brain mask from the outDir
disp(['Loading brain Mask data from ' BMfile '...']);
brainMask = readFileNifti(BMfile);
xform = brainMask.qto_xyz;
mmPerVox = brainMask.pixdim;
brainMask = logical(brainMask.data);




%% Handle the Freesurfer Segmentation

% Load the segmentation nifti file - output from freesurfer
fs = readFileNifti(freesurfer);

% Dimenstionality check 
if fs.dim(1)==size(T1,1) && fs.dim(2)==size(T1,2) && fs.dim(3)==size(T1,3)
else
    bb = mrAnatXformCoords(xform,[1 1 1;size(T1)]);
    fs.data = mrAnatResliceSpm(double(fs.data),inv(fs.qto_xyz),bb,mmPerVox,1);
end

fs = double(fs.data);

% Do some more dimensionality checks between the T1 and the freesurfer
% segmentation
if size(fs,1)==size(T1,1)+1;
    disp('freesurfer is different in size from the T1 we will clip the extra x voxels and hope its right')
    fs1 = fs; clear fs
    fs(1:size(T1,1),:,:)=fs1(1:size(T1,1),:,:); clear fs1;
end

if size(fs,2)==size(T1,2)+1;
    disp('freesurfer is differe in size from the data we clip the extra y voxels. we hope it right')
    fs1=fs;clear fs
    fs(:,1:size(T1,2),:)=fs1(:,1:size(T1,2),:); clear fs1;
end

if size(fs,3)==size(T1,3)+1;
    disp('freesurfer is differe in size from the data we clip the extra z voxels. we hope it right')
    fs1=fs;clear fs
    fs(:,:,1:size(T1,3))=fs1(:,:,1:size(T1,3)); clear fs1;
end

if size(fs,1)~=size(T1,1) || size(fs,2)~=size(T1,2) || size(fs,3)~=size(T1,3)
    error('The freesurfer segmentation file does not match the T1 data size!')
end


%% Calculate the masks using the FS labels

CSF = zeros(size(brainMask));

% ROI={'Right-Lateral-Ventricle' 'Left-Lateral-Ventricle' 'Left-Inf-Lat-Vent' 'Right-Inf-Lat-Vent' 'CSF'}
% %%http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT

CSF(fs==43) = 1; % Right-Lateral-Ventricle
CSF(fs==4)  = 2; % Left-Lateral-Ventricle
CSF(fs==5)  = 3; % Left-Inf-Lat-Vent
CSF(fs==44) = 4; % Right-Inf-Lat-Vent
CSF(fs==24) = 5; % CSF
CSF(~brainMask) = 0;

% CSF1=CSF & T1>3.8 & T1< 5;

% Create a CSF Mask
CSF = CSF & T1>4. & T1< 5;
CSFfileFS = fullfile(outDir, 'csf_seg_T1.nii.gz');
dtiWriteNiftiWrapper(single(CSF), xform, CSFfileFS);

mask = zeros(size(brainMask));
mask = double(mask);
mask(find(CSF)) = 1;
 
% Set up the white matter mask     
wm = zeros(size(brainMask));
wm(fs==2)  = 1;
wm(fs==41) = 1;
wm         = logical(wm);
[d dd] = ksdensity(T1(wm),(min(T1(wm)):0.01:max(T1(wm))));

M = dd(find(d==max(d)));
% M = mean(T1(wm)) ;
% S = std(T1(wm));

wm = wm &  T1>(M-0.03) & T1<(M+0.03);
mask(find(wm)) = 2;

% Create Tissue Mask from those FS label values above 1000
cortex = fs>1000;
[d dd] = ksdensity(T1(cortex), (min(T1(cortex)):0.01:max(T1(cortex))) );
M      = dd(find(d==max(d))); 
cortex = cortex & T1 > (M-0.03) & T1 < (M+0.03);

mask(find(cortex)) = 3;
fileFS = fullfile(outDir,'T1w_tissue.nii.gz');
dtiWriteNiftiWrapper(single(mask), xform, fileFS);

 
mask = zeros(size(brainMask));
mask(fs==43) = 1; % Right-Lateral-Ventricle
mask(fs==4)  = 1; % Left-Lateral-Ventricle
mask(fs==5)  = 1; % Left-Inf-Lat-Vent
mask(fs==44) = 1; % Right-Inf-Lat-Vent
mask(fs==24) = 1; % CSF
mask(fs==2)  = 1;
mask(fs==41) = 1;
mask(fs==2)  = 3; %WM
mask(fs==41) = 3; %WM
mask(fs>1000) = 2; %GM
segfile=fullfile(outDir,'t1_bet_seg.nii.gz');
dtiWriteNiftiWrapper(single(mask), xform, segfile);


return
 






