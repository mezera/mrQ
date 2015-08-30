function opt=mrQ_PD2WF_step5(opt,csffile,segfile)

% opt=mrQ_PD2WF_step5(opt,csffile,segfile)
%
% This is Step 6 of 6 (including Step 0) in the pipeline to build the WF
% (water fraction) map. In this step, the WF map is constructed from the
% uncalibrated PD map using the segmentation images "csffile" and
% "segfile".
%
%  ~INPUTS~
%              opt:   The "opt" structure of optimized parameters
%          csffile:   The CSF segmentation file. If no csffile is selected,
%                              the default is to take from opt's outDir, 
%                              csf_seg_T1.nii.gz
%          segfile:   The T1-weighted tissue segmentation file. If no
%                              segfile is selected, the default is to take 
%                              from opt's outDir,T1w_tissue.nii.gz            
%
%  ~OUTPUTS~
%              opt:   The updated opt structure of optimized parameters
%
% See also: mrQ_buildPD_ver2
%           Step_0: none
%           Step_1: mrQ_CalBoxPD_step1a
%           Step_2: mrQ_ScaleBoxes_step2
%           Step_3: mrQ_BoxJoinBox
%           Step_4: mrQ_smoothGain_step4b
%
% AM (C) Stanford University, VISTA

%% I. Load files 
% Get the CSF (ventricles) full-tissue segmetation and get the PD maps to
% calcultate the water fraction

if notDefined('csffile'); csffile = fullfile(opt.outDir, 'csf_seg_T1.nii.gz');end
if notDefined('segfile'); segfile = fullfile(opt.outDir, 'T1w_tissue.nii.gz');end

if (exist(csffile,'file'));
    fprintf(['Loading CSF data from: ' csffile '\n']); 
    CSF = readFileNifti(csffile);
    CSF=double(CSF.data);
else
    error(['Error, cannot find the file: '  csffile]);
end

if(exist(segfile,'file'))
    fprintf(['Loading segmentation data from: ' segfile '\n'])
    seg = readFileNifti(segfile);
    seg=double(seg.data);
else
    error(['Error, cannot find the file: '  segfile])
end

PD=readFileNifti(opt.PDfile);
xform=PD.qto_xyz;
PD=PD.data;

%% II. Calculate the CSF PD, and calibrate

% Find the white matter mean PD value from segmentation.

% Make sure the CSF ROI has reasonable PD values.  

% The CSF ROI is a result of a segmentation algorthm run on the T1-weighted
% image and cross section with T1 values. Yet the ROI may have some
% contaminations or segmentation faults, so we will create some low and
% high boundaries. Thus, PD values that are equivalent to that of the white
% matter (too low) or are double that of the white matter (too high) will
% not be considered CSF.
wmV=mean(PD(seg==3 & PD>0)); 

CSF1=CSF & PD>wmV & PD< wmV*2;

%% III. Calibrate the PD by the PD of the CSF ROI
% To calibrate the PD, we find the scalar that shifts the CSF ROI values to
% be equal to 1, that is, PD(CSF)=1.  To find the scalar, we look at the
% histogram of PD values in the CSF, focusing on its maximum. Since it's
% not trivial to find the peak, we compute the kernel density (or
% distribution estimates). For additional information, see ksdensity.m
% (standard Matlab function).

[csfValues csfDensity]= ksdensity(PD(find(CSF1)), [min(PD(find(CSF1))):0.001:max(PD(find(CSF1)))] );
CalibrationVal= csfDensity(find(csfValues==max(csfValues))); %median(PD(find(CSF)));

WF=PD./CalibrationVal(1);

% Removing outliers
WF(WF<0)=0; %too low
WF(WF>2)=2; %too high

%% IV. Save the WF map
WFfile=fullfile(opt.outDir,'WF_map.nii.gz');
dtiWriteNiftiWrapper(single(WF), xform, WFfile);
opt.WFfile=WFfile;
   save(opt.logname,'opt')
