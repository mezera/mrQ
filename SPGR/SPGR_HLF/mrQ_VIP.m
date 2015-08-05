function [AnalysisInfo]=mrQ_VIP(mrQ,outDir,WFfile,T1file,mField,T1freeval,Fullerton)
% 
% [AnalysisInfo]= mrQ_VIP(mrQ,outDir,WFfile,T1file,mField,T1freeval,Fullerton)
%
% This function loads the T1 and WF maps and calculates VIP and SIR maps.
%
% This function uses the model developed by A.M. by estimating the T1 map
% in 0.5T, 1.5T and 3T scanners for 3 subjects.
%
% This function makes the transformation from T1 and WF maps to maps of VIP
% (volume of water-interacting protons) and SIR (the water surface
% interaction ratio). The transformation is based on the fast-exchange
% model of T1 and is described in Mezer et al. (2012), in which we
% calculated the VIP, SIR and TV. The model was fit with T1 in 1.5T and 3T
% and was tested also at 0.5T (this won't work for different magnetic field
% data).
%
% 
% ABOUT THE MODEL:
% 
%   The T1 value is modeled as a weighted sum of two fast exchanging pools;
%   a free pool (with T1f = ~4.3 sec) and a hydration pool (T1h).
%       
%       1/T1 = fh/T1h + (1-fh)/T1f
% 
%   T1h is estimated as a linear function of the magnetic field (Fullerton,
%   1984). (calculated in-vitro with different tissue types)
% 
%       T1h = 1.83 x f + 25, 
% 
%   OR our new model: (calculated with in-vivo brain data)
% 
%       T1h = 0.934 x f + 93.03
% 
%          where f is the Larmor frequency for the given magnetic field. 
%   * Model estimation std dev for our values are 0.0252 and 1.8035
% 
%   Rearranging the equation above, the water fraction (fh) is given by: 
%       
%       fh = (1/T1-1/T1f) x (1/T1h-1/T1f).
% 
% 
% INPUTS:
%       mrQ       - The mrQ structure
%
%       outDir    - Directory containing the aligned SPGR data.
% 
%       WFfile    - The Water Fraction file from mrQ_WF
% 
%       T1file    - The T1 fit NIfTI
% 
%       mField    - The strength of the magnetic field. Default taken from 
%                   the dicom header. 
% 
%       T1freeval - The T1 value for free water [Default is 4.3 seconds]
% 
%       Fullerton - Boolean: 1 = use the Fullerton model for the
%                   calculation; 0 = don't use it. 
% 
% 
% OUTPUTS:
%   AnalysisInfo  - an information structure
%
%   maps that are saved
% if Fullerton==1
%     dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fitFullerton.nii.gz'));
%     dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_fitFullerton.nii.gz'));
%     
% else
%     
%    dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fit.nii.gz'));
%     dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_map.nii.gz'));
%     dtiWriteNiftiWrapper(single(TV), xform, fullfile(outDir,'TV_map.nii.gz'));
%     dtiWriteNiftiWrapper(single(SIR), xform, fullfile(outDir,'SIR_map.nii.gz'));
%
% 
% (C) Stanford University, VISTA Lab
% 


%% I. CHECK INPUTS

if notDefined('outDir')
    outDir=mrQ.spgr_initDir;
end

if(exist('Fullerton','var') && Fullerton == 1)
    disp('Using Fullerton model');
    % Using Fullerton parameters from Fullerton et al., 1984
else
    Fullerton = 0;
end

% Set the free value of T1
if(~exist('T1freeval','var') || isempty(T1freeval))
    disp('The water T1 used is 4.3 seconds ');
    % This comes from the literature (MRM 1984?)
    T1freeval = 4.3;
end


% Get the T1 file
if(exist('T1file','var') && ~isempty(T1file))
    disp(['Loading T1 data from ' T1file '...']);
else
   [ T1file,~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
end

%         T1file1= fullfile(outDir,'maps/T1_map_lsq.nii.gz');
%     T1file= fullfile(outDir,'T1_map_lsq.nii.gz');
%     if(exist(T1file,'file'))
%         disp(['Loading T1 data from ' T1file '...']);
%    
%     elseif(exist(T1file1,'file') &&  ~exist(T1file,'file') )
%           disp(['Loading T1 data from ' T1file1 '...']);
%           T1file=T1file1;
%         
%     else
%         T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',outDir);
%         if isempty(T1file)
%             error('User cancelled.')
%         else
%             disp(['Loading T1 data from ' T1file '...']);
%         end
%     end
% end


% Get the Water Fraction File
if(exist('WFfile','var') &&  ~isempty(WFfile))
    disp(['Loading WF data from ' WFfile '...']);
else
    WFfile = fullfile(outDir,'WF_map.nii.gz');
    WFfile1 = fullfile(outDir,'maps/WF_map.nii.gz');
    
    disp(['Trying  to load WF from ' WFfile '...']);
    
    if(exist(WFfile,'file'))
        disp(['Loading WF data from ' WFfile '...']);
    elseif(exist(WFfile1,'file') && ~exist(WFfile,'file'))
        disp(['Loading WF data from ' WFfile1 '...']);
        WFfile=WFfile1;
    else
        WFfile = mrvSelectFile('r','*.nii.gz','Select WF file',outDir);
        if isempty(WFfile)
            error('User cancelled.')
        else
            disp(['Loading WF data from ' WFfile '...']);
        end
    end
    tmp = 0;
end


%% II. LOAD DATA

% Load the T1
T1    = readFileNifti(T1file);
xform = T1.qto_xyz;
T1    = double(T1.data);

% Load the Water Fraction NIfTI
WF  = readFileNifti(WFfile);
WF  = double(WF.data);
WF(WF>1)=1;

infofile = fullfile(outDir,'AnalysisInfo.mat');
load(infofile);

% Set field strength
if (~exist('mField','var')|| isempty(mField))
    mField = mrQ.fieldstrength;
end

% Save files used and params for future reference
AnalysisInfo.T1forVIP  = T1file;
AnalysisInfo.WFforVIP  = WFfile;
AnalysisInfo.T1freeval = T1freeval;
AnalysisInfo.VIPdate   = date;
AnalysisInfo.fieldStrength_forVIP = mField;

save(infofile,'AnalysisInfo');


%% III. SWITCH on magnet field strength and set the larmour frequency (L)

switch mField    
    case 3
        L = 127.74;      
    case 1.5
        L = 63.87;
    case .5
        L = 21.29;       
end

if notDefined('L')
    display('Unknown magnetic field.')
    return
end


%% IV. VIP Calculation
%% IV-a. Choose whether or not to use the Fullerton model

disp('Calculating VIP ... ');

mask = (find(WF));
fh   = zeros(size(T1));

if Fullerton == 1 
    % Fullerton model
    fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(1.83.*L + 25.02)-1/T1freeval);
else
  % fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(0.934.*L+ 93.38)-1/T1freeval);
    fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(0.934.*L+ 93.03)-1/T1freeval); 
    % THIS is the right one. For some reason, I used the other one, so
    % better to redo it.
end
% fh is the fraction of interacting protons. 
% T1 is the weighted sum of the proton population under the fast exchange system.
% 1/T1= fh*(1/T1h)+ (1-fh)1/T1free.

%% IV-b. Threshold 
% Set values that are too high, too low, Inf, or Nan all to be zero.
fh(isnan(fh))   = 0;
fh(isinf(fh))   = 0;
fh(find(fh>.8)) = 0;
fh(find(fh<0))  = 0;

%% IV-c. The VIP map is calculated from VIP= fh x WF.
% VIP makes fh be a fraction of the voxel (volume).
VIP       = zeros(size(T1));
VIP(mask) = WF(mask).*fh(mask);

%% IV-d. Clean noise in the VIP signal
% Set values that are too high, too low, Inf, or Nan all to be zero.
VIP(isnan(VIP))   = 0;
VIP(isinf(VIP))   = 0;
VIP(find(VIP>.5)) = 0;
VIP(find(VIP<0))  = 0;

%% V. TV and SIR
% Calculate TV
TV=1-WF;
% Calculate the volume ratio (SIR)
SIR=VIP./TV;

%% VI. Save Output

if Fullerton==1
    dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fitFullerton.nii.gz'));
    dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_fitFullerton.nii.gz'));
    dtiWriteNiftiWrapper(single(TV), xform, fullfile(outDir,'TV_map.nii.gz'));
    
else
    
    dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fit.nii.gz'));
    dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_map.nii.gz'));
    dtiWriteNiftiWrapper(single(TV), xform, fullfile(outDir,'TV_map.nii.gz'));
    dtiWriteNiftiWrapper(single(SIR), xform, fullfile(outDir,'SIR_map.nii.gz'));
    
end

return