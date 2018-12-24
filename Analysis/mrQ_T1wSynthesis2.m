function [saveName,saveName1,saveName2,saveName3] =mrQ_T1wSynthesis2(mrQ,WFfile,T1file,BMfile, ...
                         outDir,symTR,symFA, saveName,saveName1,saveName2,saveName3,FullBMfile)
% [saveName,saveName1] =mrQ_T1wSynthesis1(mrQ,WFfile,T1file,BMfile, ...
%                        outDir,symTR,symFA, saveName,saveName1,saveName2,saveName3,FullBMfile)
%
% This function creates a series of synthetic T1-weighted images and saves
% them to the dataDir. One set of T1w images accounts for the PD, taken
% from the water fraction file, while the other set of T1w images does not.
% 
% ~INPUTS~
%           mrQ:   The mrQ structure
%        WFfile:   The path to the directory where the water fraction map
%                      is located. If unavailable, the default will take
%                      the linearly fitted M0 map.
%        T1file:   The path to the directory where the align.mat file 
%                      exists from the getSEIR function earlier.
%        BMFile:   The path to the directory where the brain mask is
%                      located. It will be created if it doesn't yet exist.
%        outDir:   The path to where you would like to save the data.
%         symTR:   Simulated TR. [Default is 30 milliseconds]
%         symFA:   Simulated flipAngles. [Default is 30 degrees]
%      saveName:   The path to the directory where the synthetic T1
%                      weighted image will be saved.
%     saveName1:   The path to the directory where the synthetic T1
%                      weighted image (without accounting for PD) will be
%                      saved.
%    FullBMfile:   The path to the directory where the full brain mask is
%                      located. It will be created if it doesn't yet exist.
%
% ~OUTPUTS~
%      saveName:   The path to the directory where the synthetic T1
%                      weighted NIfTI image was saved.
%     saveName1:   The path to the directory where the synthetic T1
%                      weighted NIfTI image (without accounting for PD) was
%                      saved.
%    saveName2:   The path to the directory where the synthetic T1
%                      weighted NIfTI image (with accounting for 1./PD to bosththe contrast) was
%                      saved.
% See also:
%   mrQ_fitT1M0.m 
% 
% (C) Stanford University, VISTA Lab [2012]
%(C) Hebruew University, Mezer Lab [2018]
%
%

%% I. Check INPUTS and load files

if(exist('T1file','var') && ~isempty(T1file))
    disp(['Loading T1 data from ' T1file '...']);
else
    [T1file,~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
  %  T1file= fullfile(mrQ.spgr_initDir,'T1_map_lin.nii.gz');
end

t1=readFileNifti(T1file);t1=t1.data;

if(exist('WFfile','var') && ~isempty(WFfile))
    disp(['Loading PD data from ' WFfile '...']);
elseif isfield(mrQ,'maps')
    WFfile=mrQ.WFfile;
else
    % If there is no water fraction map, default will take the linearly
    % fitted M0 map, which was made in the function mrQfit_T1M0_Lin.m.
   [~, WFfile,~]=mrQ_get_T1M0_files(mrQ,0,1,0);

end

PD=readFileNifti(WFfile);PD=PD.data;

if notDefined('outDir')
    outDir=mrQ.outDir;
end

if ~exist(fullfile(outDir,'SyntheticT1w'),'dir')
    mkdir(fullfile(outDir,'SyntheticT1w'));
end

if (~exist('saveName_a','var') || isempty(saveName_a)),
saveName_a=fullfile(outDir,'SyntheticT1w','T1w_a.nii.gz');
end

if (~exist('saveName1_a','var') || isempty(saveName1_a)),
saveName1_a=fullfile(outDir,'SyntheticT1w','T1w1_a.nii.gz');
end

if (~exist('saveName2_a','var') || isempty(saveName2_a)),
saveName2_a=fullfile(outDir,'SyntheticT1w','T1w2_a.nii.gz');
end

if (~exist('saveName3_a','var') || isempty(saveName3_a)),
saveName3_a=fullfile(outDir,'SyntheticT1w','T1w3.nii.gz');
end

if (~exist('saveName3_b','var') || isempty(saveName3_b)),
saveName3_b=fullfile(outDir,'SyntheticT1w','T1w4.nii.gz');
end


if (~exist('saveName','var') || isempty(saveName)),
saveName=fullfile(outDir,'SyntheticT1w','T1w.nii.gz');
end

if (~exist('saveName1','var') || isempty(saveName1)),
saveName1=fullfile(outDir,'SyntheticT1w','T1w1.nii.gz');
end

if (~exist('saveName2','var') || isempty(saveName2)),
saveName2=fullfile(outDir,'SyntheticT1w','T1w2.nii.gz');
end

if (~exist('saveName3','var') || isempty(saveName3)),
saveName3=fullfile(outDir,'SyntheticT1w','T1w3.nii.gz');
end





if (~exist('symTR','var') || isempty(symTR)),
    symTR = 30;
end

if (~exist('flipAngleIn','var') || isempty(symFA)),
    symFA = 30;
end

if(exist('BMfile','var') && ~isempty(BMfile))
     disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    xform=brainMask.qto_xyz;

    brainMask=logical(brainMask.data);
else
    % Look for and load the brain mask - create one if necessary
   [~,~,BMfile]=mrQ_get_T1M0_files(mrQ,0,0,1);
%    BMfile = fullfile(mrQ.spgr_initDir,'HeadMask.nii.gz');
end

if (~exist('FullBMfile','var') || isempty(FullBMfile)),
    % Look for and load the brain mask - create one if necessary
        FullBMfile=BMfile;
end

 mask=readFileNifti(FullBMfile);
 mask=mask.data;

%% II. Calculate the synthetic T1 images

t1w = zeros(size(t1));

% Calculate values for t1 and fa
t1 = t1.*1000; % msec
fa = symFA./180.*pi;

% Calculate the synthetic t1 images 
t1w = PD.*( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa)));

% ... For future exploration ...
%t1w = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa))); % no PD
% t1w = t1w.*(1-PD);  % get the PD to be in our side 

% scale up to have big number like a typical MRI
% t1w = t1w.*10000;

% % clip outliers
% M=mean(t1w(brainMask));
% S=std(t1w(brainMask));
% 
% up=min(10000,M+3*S);
% down=max(0,M-3*S);
%  t1w(t1w<down)=down;
%  t1w(t1w>up)=up;    
% 
%   t1w(isnan(t1w))=down;
%  
%  t1w(isinf(t1w))=up;
 % clip what we define as notbrain
 
 fprintf('\n Calculating T1w...              \n');

  t1w(~mask)=0;
  
%% III. Now calculate the synthetic T1 images, without PD
  t1ww = zeros(size(t1));

% Calculate values for t1 and fa
fa = symFA./180.*pi;

% Calculate the synthetic t1 images 
t1ww = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa)));
%    ^^ not multiplying by PD here!


% for future exploration
%t1w = ( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa))); % no PD
% t1w = t1w.*(1-PD);  % get the PD to be in our side 

% scale up to have big number like a typical MRI
% t1ww = t1ww.*10000;

% clip outlayers
% M=mean(t1ww(brainMask));
% S=std(t1ww(brainMask));
% 
% up=min(10000,M+3*S);
% down=max(0,M-3*S);
%  t1ww(t1ww<down)=down;
%  t1ww(t1ww>up)=up;    
% 
%    t1ww(isnan(t1ww))=down;
%  
%  t1ww(isinf(t1ww))=up;
 
 % clip what we define as notbrain
  t1ww(~mask)=0;
  
  
%% Iv.  Lets overcome noise values using PD
RatioM=median(t1w(logical(mask))./t1ww(logical(mask)));
Ratio=(t1w./t1ww.*RatioM);
t1ww1=t1ww;
t1ww1(Ratio<0.2)=t1w(Ratio<0.2)/RatioM;
t1ww1(Ratio>1.9)=t1w(Ratio>1.9)/RatioM;

%% v. Now calculate the synthetic T1 images, with 1./PD
  t1www = zeros(size(t1));

% Calculate values for t1 and fa
fa = symFA./180.*pi;

% Mmake PD contrast to be similar to T1w image (1/PD).
PD1=1./PD;
%find values that are too small
% MinVal=prctile(t1ww(logical(mask)),10);
% MinVal1=prctile(PD1(logical(mask)),1);
% PD1(isinf(PD1))=1;
% PD1(PD1>=2.5 | PD1<=0.02 |t1<50)=PD(PD1>=2.5 | PD1<=0.02 | t1<50); % 
t1www=PD1.*t1ww;
  t1www(~mask)=0;

RatioM=median(t1w(logical(mask))./t1www(logical(mask)));
Ratio=t1w./(t1www.*RatioM);
t1www1=t1www;
t1www1(Ratio<0.1)=t1w(Ratio<0.1)./RatioM;
t1www1(Ratio>1.9)=t1w(Ratio>1.9)./RatioM;
    
  
% % Calculate the synthetic t1 images 
% t1www = 1./PD.*( (1-exp(-symTR./t1)).*sin(fa)./(1-exp(-symTR./t1).*cos(fa)));
% %    ^^  multiplying by 1./PD here will bost the contrat!
% MinVal1=prctile(T1WRaw1(logical(mask)),10);
% MinVal1=prctile(T1WRaw1(logical(mask)),10);
% 
% % this image have problem with dividing small values by small values. this can solve some of it.
% %PD1(PD1<0.3)=1;
% %
%  
%  % clip what we define as notbrain
%   t1www(~mask)=0;
% 
%   
% %%

outFile  = fullfile(mrQ.InitSPGR.spgr_initDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);
    flipAngles = [s(:).flipAngle];
Hfa=find(flipAngles>10);
Lfa=find(flipAngles<10);

T1WRaw=mean(cat(4,s(Hfa).imData),4);
MinVal=prctile(T1WRaw(logical(mask)),10);
T1WRaw1=mean(cat(4,s(Lfa).imData),4);
MinVal1=prctile(T1WRaw1(logical(mask)),10);

T1WRaw1(T1WRaw1<=MinVal1 & T1WRaw<=MinVal)=100;

t1wwww=T1WRaw./T1WRaw1;
  

  RatioM=median(T1WRaw(logical(mask))./t1wwww(logical(mask)));
Ratio=T1WRaw./(t1wwww.*RatioM);
t1wwww1=t1wwww;
t1wwww1(Ratio<0.1)=T1WRaw(Ratio<0.1)./RatioM;
t1wwww1(Ratio>1.9)=T1WRaw(Ratio>1.9)./RatioM;
  

%% VI. Save out the resulting NIfTI files

% Write them to disk
disp('2.  Saving synthetic T1w images.');
dtiWriteNiftiWrapper(single(t1w), mrQ.xform, saveName);
dtiWriteNiftiWrapper(single(t1ww), mrQ.xform, saveName1);
dtiWriteNiftiWrapper(single(t1www), mrQ.xform, saveName2);
dtiWriteNiftiWrapper(single(t1wwww), mrQ.xform, saveName3);

%more styT1 options
% dtiWriteNiftiWrapper(single(t1w1), mrQ.xform, saveName_a);
% dtiWriteNiftiWrapper(single(t1ww1), mrQ.xform, saveName1_a);
% dtiWriteNiftiWrapper(single(t1www1), mrQ.xform, saveName2_a);
% dtiWriteNiftiWrapper(single(t1wwww1), mrQ.xform, saveName3_a);
% dtiWriteNiftiWrapper(single(T1WRaw), mrQ.xform, saveName3_b);


return