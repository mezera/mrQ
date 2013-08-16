function [AnalysisInfo]=mrQ_VIP(outDir,WFfile,T1file,mField,T1freeval,Fullerton)
% 
% mrQ_VIP(outDir,WFfile,T1file,mField,T1freeval,Fullerton)
%# Load T1 and WF maps and calculate the VIP and SIR map using the model develop by
% A.M BY ESTIMATING t1 MAP IN 0.5T 1.5T & 3T FOR 3 SUBJECTS
%
% in this function we will make the transformation from T1 and Water fraction map to VIP
% voulume of water interacting proton.
%the transformation is based on fast exgance model of T1 and describe in
%mezer et.al. 2012. the model was fit with T1 in 1.5T and 3T) and was
%testesd also to 0.5T this won't work for different magnetic feild data.
%we also calclate the SIR the water surface interaction ratio. we calculate
%the VIP/TV see the same article.
%
 
% ABOUT THE MODEL:
% 
%   The T1 value is modeled as a weighted sum of two fast exchanging pools; a
%   free pool (with T1f = ~4.3 sec) and a hydration pool (T1h).
%       
%       1/T1 = fh/T1h + (1-fh)/T1f
% 
%   T1h is estimated as a linear function of the magnetic field (Fullerton
%   1984). (calculate in-vitro with different tissue types)
% 
%       T1h = 1.83 x f + 25, 
% 
%   OR our new model: (calculate with in-vivo brain data)
% 
%       T1h = 0.934 x f + 93.03
% 
%   Where f is the Larmor frequency for the given magnetic field. 
%   * Model estimation std for our values are 0.0252 and 1.8035
% 
%   Rearranging the equation above, the water fraction (fh) is given by: 
%       
%       fh = (1/T1-1/T1f) x (1/T1h-1/T1f).
% 
% 
% INPUTS:
%       outDir    - Directory containing the aligned SPGR data.
% 
%       WFfile    - The Water Fraction file from mrQ_WF
% 
%       T1file    - The T1 fit nifti
% 
%       mField    - The strength of the magnetic field. Default taken from 
%                   the dicom header. 
% 
%       T1freeval - Defaults to 4.3
% 
%       Fullerton - Boolean: 1 = use the Fullerton model for the
%                   calculation, 0 = don't use it. 
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
% (C) Stanford University, VISTA Lab
% 


%% CHECK INPUTS

if(exist('Fullerton','var') && Fullerton == 1)
    disp('Using Fullerton model');
    % Using Fullerton parameters form Fullerton et.al. 1984
else
    Fullerton = 0;
end

% Set the free value of T1
if(~exist('T1freeval','var') || isempty(T1freeval))
    disp('The water T1 used is 4.3sec ');
    % this comes from the litrature MRM 1984?
    T1freeval = 4.3;
end


% Get the T1 file
if(exist('T1file','var') && ~isempty(T1file))
    disp(['Loading T1 data from ' T1file '...']);
else
    T1file= fullfile(outDir,'T1_map_lsq.nii.gz');
    if(exist(T1file,'file'))
        disp(['Loading T1 data from ' T1file '...']);
    else
        T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',outDir);
        if isempty(T1file)
            error('User cancelled.')
        else
            disp(['Loading T1 data from ' T1file '...']);
        end
    end
end


% Get the Water Fraction File
if(exist('WFfile','var') &&  ~isempty(WFfile))
    disp(['Loading WF data from ' WFfile '...']);
else
    WFfile = fullfile(outDir,'WF_map.nii.gz');
    disp(['trying  to load WF from ' WFfile '...']);
    if(exist(WFfile,'file'))
        disp(['Loading WF data from ' WFfile '...']);
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



%% LOAD DATA

% Load the T1
T1    = readFileNifti(T1file);
xform = T1.qto_xyz;
T1    = double(T1.data);

% Load the Water Fraction Nifti
WF  = readFileNifti(WFfile);
WF  = double(WF.data);
WF(WF>1)=1;


infofile = fullfile(outDir,'AnalysisInfo.mat');
load(infofile);

% Set field strength
if (~exist('mField','var')|| isempty(mField))
    mField = AnalysisInfo.fieldStrength;
end

% Save files used and params for future reference
AnalysisInfo.T1forVIP  = T1file;
AnalysisInfo.WFforVIP  = WFfile;
AnalysisInfo.T1freeval = T1freeval;
AnalysisInfo.VIPdate   = date;
AnalysisInfo.fieldStrength_forVIP = mField;

save(infofile,'AnalysisInfo');


%% SWITCH on magnet field strength and set the larmour frequency (L)

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


%% VIP Calculation

disp('Calculating VIP ... ');

mask = (find(WF));
fh   = zeros(size(T1));

if Fullerton == 1 
    % Fullerton model
    fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(1.83.*L + 25.02)-1/T1freeval);
else
    % fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(0.934.*L+ 93.38)-1/T1freeval);
    
    fh(mask) = (1./T1(mask)- 1/T1freeval)./(1000./(0.934.*L + 93.03)-1/T1freeval); % THIS is the right one for some reson i used the other so better to redo it
end


%% Threshold 

fh(isnan(fh))   = 0;
fh(isinf(fh))   = 0;
fh(find(fh>.8)) = 0;
fh(find(fh<0))  = 0;


%% The VIP map is calculated from VIP= fh x WF.

VIP       = zeros(size(T1));
VIP(mask) = WF(mask).*fh(mask);


%% Clean noise

VIP(isnan(VIP))   = 0;
VIP(isinf(VIP))   = 0;
VIP(find(VIP>.5)) = 0;
VIP(find(VIP<0))  = 0;



%% TV and SIR
%calculate TV
TV=1-WF;
%calculate the voulume ratio SIR
SIR=VIP./TV;



%% Save Output

if Fullerton==1
    dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fitFullerton.nii.gz'));
    dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_fitFullerton.nii.gz'));
    
else
    
    dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,'T1wVIP_fit.nii.gz'));
    dtiWriteNiftiWrapper(single(VIP), xform, fullfile(outDir,'VIP_map.nii.gz'));
    dtiWriteNiftiWrapper(single(TV), xform, fullfile(outDir,'TV_map.nii.gz'));
    dtiWriteNiftiWrapper(single(SIR), xform, fullfile(outDir,'SIR_map.nii.gz'));

end

return