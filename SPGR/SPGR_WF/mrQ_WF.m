function [mrQ,WFfile,CalibrationVal] = mrQ_WF(mrQ,dataDir,T1file,PDfile,Gainfile,B1file)
% [mrQ,WFfile,CalibrationVal] = mrQ_WF(mrQ,dataDir,T1file,PDfile,Gainfile,B1file);
% 
% 
% This function calculates the constant deviation of the PD. Since some
% flipangles have higher noise than others, it finds the flipangles with
% the highest signal of the rawaligned data within the CSF, (as the highest
% signal should reflect the highest SNR). Than it creates an ROI which is
% the CSF only within a limited ara around the ventricles, while excluding
% voxel with outliers in PD or in the raw data.
% Next PD is calculated in the ROI, annd assuming the value of it should
% be 1,  the correction constant of the PD, and apply that constant to the
% PD image to create a WF map.
% 

% Input:
%         mrQ 
%         dataDir -   where the structure with the raw data is located.
%                     default is mrQ.spgr_initDir 
%         T1file -    the path to the T1 map. This is used to create the 
%                     CSF ROI. default is using mrQ_get_T1M0_files 
%         PDfile -    the path to the PD map.  default is: opt.PDfile 
%         Gainfile -  the path to a Gain map. This is used for the PD
%                     calculation. the gain is defined differently,
%                     depending on the M0 calculation (if it's a multiCoil
%                     fit than we take a different map) 
%         B1file -    the path to a B1 map, this will be used to correct
%                     the FA value. default is mrQ.B1FileName
%% I : check input
if notDefined('dataDir');
    dataDir = mrQ.spgr_initDir;
end

if notDefined('T1file')
    [ T1file,~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
end
T1=readFileNifti(T1file); T1=T1.data;

if notDefined('PDfile')
    load(mrQ.opt_logname);
    PDfile=opt.PDfile;
end
PD=readFileNifti(PDfile);PD=PD.data;

if notDefined('B1file')
    B1file=mrQ.B1FileName;
end
 B1=readFileNifti(B1file);B1=double(B1.data);
 
 % Gain file will be loaded andchecked later. When M0 was calculated using
 % the multicoil information we'd like to only load the Gain file matching
 % to our flip angle of choice. 

% currently, flipangle is chosen based on signal witin CSF. 
%  if instead it would be decided based on signal with WM, that we would
%  need a brainmask: 

% if notDefined('BMfile')
%     BMfile= mrQ.BrainMask;
% end
% BM=readFileNifti(BMfile); BM=BM.data;

%% II. Load aligned raw data

% next we load the file that contains the aligned raw data
% this file contains: xform, mmPerVox, and the struct 's'
% s.flipAngle contains the flipngle, and s.imData the data.
outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data
disp(['Loading aligned data from ' outFile '...']);
load(outFile);

%% III. find area of high T1=[4 - 5]  around the vetricles

WFmask=ones(size(T1));

% clip around the ventricals:
if notDefined('boxsize')
    boxsize(1)=30;
    boxsize(2)=40;
    boxsize(3)=20;
end
sz=size(WFmask); szH=round(sz./2);
XX=boxsize(1)./round(mmPerVox(1));
YY=boxsize(2)./round(mmPerVox(2));
ZZ=boxsize(3)./round(mmPerVox(3));

WFmask(szH(1)+XX:end,:,:)=0;
WFmask(1:szH(1)-XX,:,:)=0;

WFmask(:,1:szH(2)-YY,:)=0;
WFmask(:,szH(2)+YY:end,:,:)=0;

WFmask(:,:,1:szH(3)-ZZ)=0;
WFmask(:,:,szH(3)+ZZ:end)=0;

% find areas within the ventricles area with high value of T1
WFmask= WFmask & T1<=5 & T1>=4.2 ;
%% IV. data for WF:
% % loop over raw (aligned and combined) data and find the flip angle(s) with
% % the best SNR. it is hard to measure noise. so we will use the signal in WM
% % as a reference because it's a place with good signal across flip angles)

for ii=1:length(s)
    CSF_Sig(ii)=median(s(ii).imData(WFmask)); 
end

[~,SelectedFA]=(max(CSF_Sig));

%%

%  Instead of finding the best FA within CSF, another option is using WM
%  ROI, 
% %     WMmask=?
% WMmask=ones(size(T1));
% WMmask(T1<0.7)=0;
% WMmask(T1>1.1)=0;
% WMmask=WMmask & BM;
% 
% for ii=1:length(s)
%     
%     rawDat=s(ii).imData;
% 
%     %     per file??
%     %     mean?
%     
%      WM_Sig(ii)=median(s(ii).imData(WMmask)); 

% end
% 
% FAloc=find(WMdata==max(WMdata));



%% V. data for WF:

datForPD=s(SelectedFA).imData;
TR=s(SelectedFA).TR;
fa=s(SelectedFA).flipAngle;

%% VI : exclude M0 outliers in the ROI

WFmask=WFmask & PD<prctile(PD(WFmask),99) & PD>prctile(PD(WFmask),1) & datForPD< prctile(datForPD(WFmask),99) &  datForPD>prctile(datForPD(WFmask),1);


%% VII. calculate PD
% calculate  PD within the selected ROI and flip angle(s) (using the SPGR
% equation, T1, and gain 


% [note gain for multi coil is define for the sum of multi coils, and for
% the combine as arrived from the scanner])

if notDefined('Gainfile')
    if mrQ.PDfit_Method==1
        %Gain is the default ??
        load(mrQ.opt_logname);
        Gainfile=opt.Gainfile;
        
    elseif mrQ.PDfit_Method==2 || mrQ.PDfit_Method==3
        
        % multicoil M0 was computed
        % Gain is defined as the sum of multicoils.        
        Gainfile=mrQ.MultiCoilSummedFiles.(gainField);
        
    end
end
Gain=readFileNifti(Gainfile); Gain=Gain.data;

% calculate PD for the FA with the best Signal

% change FA into radians
fa = (fa.*B1)./180.*pi;

% calculate PD assuming we're in CSF (so T1 is 4300ms)
PDcsf=datForPD./(Gain.* ( (1-exp(-TR./4300)).*sin(fa)./(1-exp(-TR./4300).*cos(fa))));
 
%% VIII. scale this roi to have pd of 1. --> this is our global scale

[csfValues, csfDensity]= ksdensity(PDcsf((WFmask)), [min(PDcsf((WFmask))):0.001:max(PDcsf((WFmask)))] ); % figure;plot(csfDensity,csfValues)
CalibrationVal= csfDensity(csfValues==max(csfValues)); %median(PD(find(CSF)));
CalibrationVal=1./CalibrationVal;

%% IX. apply to PD images to make WFfile --> save

WF=PD.*CalibrationVal(1);

WFfile=fullfile(dataDir,'WF_map.nii.gz');
dtiWriteNiftiWrapper(WF,xform,WFfile);

mrQ.WFfile=WFfile;
mrQ.ScalePD_2_WF=CalibrationVal;

    

