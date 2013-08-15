function opt=mrQ_PD2WF_step5(opt,csffile,segfile)
% opt=mrQ_PD2WF_step5(opt,csffile,segfile)
% calculate Water fraction map from unclaibrated PD map using the segmentaion images
%
% AM (C) Stanford University, VISTA

%% get the CSF (ventricals) full tissue segmetation and PD maps to calcultate Water fraction (WF)
if notDefined('csffile');csffile = fullfile(opt.outDir, 'csf_seg_T1.nii.gz');end
if notDefined('segfile');segfile = fullfile(opt.outDir, 't1_bet_seg.nii.gz');end

if(exist(csffile,'file')); fprintf(['Loading CSF data from: ' csffile '\n']); CSF = readFileNifti(csffile);CSF=double(CSF.data);
else   error(['error , can not find the file: '  csffile])    ;end

if(exist(segfile,'file'));fprintf(['Loading segmentation data from: ' segfile '\n']);seg = readFileNifti(segfile);seg=double(seg.data);
else error(['error , can not find the file: '  segfile]);end

PD=readFileNifti(opt.PDfile); xform=PD.qto_xyz;PD=PD.data;

%% calcute the CSF PD

% find the white matter mean pd value from segmetation.
wmV=mean(PD(seg==2 & PD>0));

% assure thata the CSF ROI have pd value that are resnable.  The csf roi is a reslut of segmentation algoritim runed on the
% T1wighted image and cross section with T1 values. Yet  the ROI may have some contaminations or segmentation faules .
%Therefore, we create some low and up bonderies. No CSF with PD values that are the white matter PD value(too low) or double the white matter values (too high).
CSF1=CSF & PD>wmV & PD< wmV*2;

%To calibrate the PD we find the scaler that shift the csf ROI to be eqal to 1. --> PD(CSF)=1;
% To find the scale we look at the histogram of PD value in the CSF. Since it's not trivial to find the peak we compute the kernel density (or
% distribution estimates). for detail see ksdensity.m
%The Calibrain vhistogram of the PD values in the let find the scalre from the maxsimum of the csf values histogram
[csfValues csfDensity]= ksdensity(PD(find(CSF1)), [min(PD(find(CSF1))):0.001:max(PD(find(CSF1)))] );
CalibrationVal= csfDensity(find(csfValues==max(csfValues)));% median(PD(find(CSF)));

%% calibrate the pd by the pd of the csf roi
WF=PD./CalibrationVal;

% let cut outlayers
WF(WF<0)=0;
WF(WF>2)=2;

%% save the WM map
WFfile=fullfile(opt.outDir,'WF_map.nii.gz');
dtiWriteNiftiWrapper(single(WF), xform, WFfile);
opt.WFfile=WFfile;
 save(opt.logname,'opt')
