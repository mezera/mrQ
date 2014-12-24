function [saveName,saveName1] =mrQ_T1wSynthesis1(mrQ,WFfile,T1file,BMfile,symTR,symFA, saveName,saveName1,FullBMfile)
% 
% mrQ_T1wSynthesis1(dataDir,B1file,outDir,trIn,flipAngleIn)
% 
% Create a series of synthetic T1w images and save them to the dataDir.
% 
% INPUTS:
%   T1file: the path to the directory where the align.mat file exists from
%            the getSEIR function earlier.
% 
%   B1file:  Specify where the B1 inhomogeneity map exists. If you leave it
%            empty, it will use the B1 file from the data directory.
% 
%   outDir:  the path to where you would like to save the data.
% 
%   trIn:    The TR [defaults to 30]
% 
%   flipAngleIn: The flip angle [defaults to 30]
% 
% 
% See also:
%   mrQ_fitT1M0.m 
% 
% 
% (C) Stanford University, VISTA Lab [2012]

%% Check INPUTS


if(exist('T1file','var') && ~isempty(T1file))
    disp(['Loading T1 data from ' T1file '...']);
else
    T1file=mrQ_getT1file(mrQ);
  %  T1file= fullfile(mrQ.spgr_initDir,'T1_map_lin.nii.gz');
end
t1=readFileNifti(T1file);t1=t1.data;

if(exist('WFfile','var') && ~isempty(WFfile))
    disp(['Loading PD data from ' WFfile '...']);
else
    WFfile=mrQ.maps.WFpath;
    
end
PD=readFileNifti(WFfile);PD=PD.data;



if (~exist('saveName','var') || isempty(saveName)),
saveName=fullfile(mrQ.OutPutNiiDir,'T1w','T1w.nii.gz');
end
if (~exist('saveName1','var') || isempty(saveName1)),

saveName1=fullfile(mrQ.OutPutNiiDir,'T1w','T1w1.nii.gz');

end


if (~exist('symTR','var') || isempty(symTR)),
    symTR = 30;
end

if (~exist('flipAngleIn','var') || isempty(symFA)),
    symFA = 30;
end

if (~exist('BMfile','var') || isempty(BMfile)),
% Look for and load the brain mask - create one if necessary
BMfile = fullfile(mrQ.spgr_initDir,'HeadMask.nii.gz');
end

if(exist(BMfile,'file'))
    disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    xform=brainMask.qto_xyz;

    brainMask=logical(brainMask.data);

end

if (~exist('FullBMfile','var') || isempty(FullBMfile)),
    % Look for and load the brain mask - create one if necessary
    
    if isfield(mrQ.AnalysisInfo,'maskSynthesis')
        FullBMfile = mrQ.AnalysisInfo.maskSynthesis;
    else
        FullBMfile=BMfile;
    end
end

 mask=readFileNifti(FullBMfile);
 mask=mask.data;

%% III. Calculate the synthetic t1 images



t1w = zeros(size(t1));

% Calculate values for t1 and fa
t1 = t1.*1000; % msec
fa = symFA./180.*pi;

% Calculate the synthetic t1 images 
t1w = PD.*( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa)));

% for future exploration
%t1w = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa))); % no PD
% t1w = t1w.*(1-PD);  % get the PD to be in our side 

% scale up to have big number like a typical MRI
t1w = t1w.*10000;

% clip outlayers
M=mean(t1w(brainMask));
S=std(t1w(brainMask));

up=min(10000,M+3*S);
dwon=max(0,M-3*S);
 t1w(t1w<dwon)=dwon;
 t1w(t1w>up)=up;    

  t1w(isnan(t1w))=dwon;
 
 t1w(isinf(t1w))=up;
 % clip what we define as notbrain
  t1w(~mask)=0;
  
  %% now without PD
  t1ww = zeros(size(t1));

% Calculate values for t1 and fa
fa = symFA./180.*pi;

% Calculate the synthetic t1 images 
t1ww = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa)));

% for future exploration
%t1w = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa))); % no PD
% t1w = t1w.*(1-PD);  % get the PD to be in our side 

% scale up to have big number like a typical MRI
t1ww = t1ww.*10000;

% clip outlayers
M=mean(t1ww(brainMask));
S=std(t1ww(brainMask));

up=min(10000,M+3*S);
dwon=max(0,M-3*S);
 t1ww(t1ww<dwon)=dwon;
 t1ww(t1ww>up)=up;    

   t1ww(isnan(t1ww))=dwon;
 
 t1ww(isinf(t1ww))=up;
 
 % clip what we define as notbrain
  t1ww(~mask)=0;

  
  
  
%% IV. Save out the resulting nifti files

% Write them to disk
disp('2.  saving a syntetic T1w images');
dtiWriteNiftiWrapper(single(t1w), xform, saveName);
dtiWriteNiftiWrapper(single(t1ww), xform, saveName1);

return

