function [T1wfs_4file] = mrQ_T1wSynthesis(dataDir,B1file,outDir,trIn,flipAngleIn)
% 
% mrQ_T1wSynthesis(dataDir,B1file,outDir,trIn,flipAngleIn)
% 
% Create a series of synthetic T1w images and save them to the dataDir.
% 
% INPUTS:
%   dataDir: the path to the directory where the align.mat file exists from
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

if (~exist('dataDir','var') || isempty(dataDir)),
    dataDir = pwd;
end

if (~exist('outDir','var') || isempty(outDir)),
    outDir = dataDir;
end

if(~exist(outDir,'dir')),
    mkdir(outDir);
end

if (~exist('B1file','var') || isempty(B1file)),
    B1file = fullfile(dataDir, 'B1_Map.nii.gz');
end

if (~exist('trIn','var') || isempty(trIn)),
    trIn = 30;
end

if (~exist('flipAngleIn','var') || isempty(flipAngleIn)),
    flipAngleIn = 30;
end


%% I. Load aligned data

outFile = fullfile(dataDir,'dat_aligned.mat');
disp(['Loading aligned data from ' outFile '...']);
load(outFile);


% Look for and load the B1 map if it exists
if exist(B1file,'file')
    B1 = readFileNifti(B1file);
    B1 = double(B1.data);  
else
    disp( 'wanrning can not find B1 map. No B1 map will be used!!!');
    B1 = ones(size(s(1).imData));
end


% Look for and load the brain mask - create one if necessary
BMfile = fullfile(outDir,'brainMask.nii.gz');

if(exist(BMfile,'file'))
    disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    brainMask=logical(brainMask.data);
else
    brainMask=BrainMaskFit(s,mmPerVox,[],[],1,xform,BMfile);
end


%% II. Perform linear fit of T1 and M0

% Return the flip angles and TR from the aligned data structure (s)
flipAngles = [s(:).flipAngle];
tr = [s(:).TR];

% Make sure that all TRs are the same and set the TR
if(~all(tr == tr(1))), error('TR''s do not match!'); end

% Set the TR
tr = tr(1);

% Do the fitting: Fitting routine from RFD and NS
disp('1. Fiting lineary T1 and M0');
[t1,M0] = relaxFitT1(cat(4,s(:).imData),flipAngles,tr,B1);


%% III. Calculate the synthetic t1 images

% Create a white-matter mask using a range of t1 values within the original
% brain mask
wmmask    = brainMask & t1>0.850 & t1<1.050;

% Perform 3-D order-statistic filtering...
[wmmask1] = ordfilt3D(wmmask,6);
wmmask    = wmmask & wmmask1; clear wmmask1;

[params1,gains,rs] = fit3dpolynomialmodel(M0,wmmask==1,2);
wmmask             = logical(wmmask);

Imsz = size(brainMask);

% Construct a three-dimenstional polynomial matrix - KNK code 
[Poly,str] = constructpolynomialmatrix3d(Imsz,find(ones(Imsz)),2);
%   Returns <Poly>, a matrix of dimensions length(find(ones(Imsz))) x N
%   with polynomial basis functions evaluated at <find(ones(Imsz))> in
%   the columns.  the polynomial basis functions are evaluated
%   over the range [-1,1] which is presumed to correspond to
%   the beginning and ending element along each of the three dimensions.
%   (if a dimension has only one element, the values are all set to 1.)
%   also, return <str>, the algebraic expression that corresponds to
%   the columns of <Poly>.  'x' refers to the first matrix dimension; 'y'
%   refers to the second matrix dimension; 'z' refers to the third 
%   matrix dimension.

% Calculate the gain and proton density
Gain = reshape(Poly*params1',Imsz);
PD   = M0./Gain;

% Create the empty t1 images
% t1w_1 = zeros(size(t1));
% t1w_2 = zeros(size(t1));
% t1w_3 = zeros(size(t1));
t1w_4 = zeros(size(t1));

% Calculate values for t1 and fa
t1 = t1.*1000; % msec
fa = flipAngleIn.*B1;
fa = fa./180.*pi;

% Calculate the synthetic t1 images 
%t1w_1 = M0.*(1-exp(-trIn./t1)).*sin(fa)./(1-exp(-trIn./t1).*cos(fa));
t1w_2 = (1-exp(-trIn./t1)).*sin(fa)./(1-exp(-trIn./t1).*cos(fa));
%t1w_3 = t1w_2./PD;
t1w_4 = t1w_2.*PD;

% Scale the values for t1w_3
% t1w_3 = t1w_3.*10000;
% t1w_3(isinf(t1w_3)) = 0;
% t1w_3(t1w_3>1e5) = 0;

% Calculate the mean and standard deviation of t1w_3
%M = mean(t1w_3(brainMask));
%S = std(t1w_3(brainMask));

% Scale down values in t1w_3 greater than 2*SD above the mean by 0.1
%t1w_3(t1w_3>(M+2*S)) = t1w_3(t1w_3>(M+2*S)).*0.1;

% Calculate the mean and standard deviation of M0
 M = mean(M0(brainMask));
 S = std(M0(brainMask));
% 
% % Create morphological structuring element (STREL) and dilate/fill the M0
% % image
% se   = strel('ball',5,5);
% mask = zeros(size(brainMask));
% for i = 1:size(M0,3)
%     I = M0(:,:,i);
%     I2 = imdilate(I,se);
%     I2 = imfill(I2,'holes');
%     mask(:,:,i) = I2>(M-S*2);
% end
% mask = logical(mask);

% Set values outside the calculated mask to zero
%t1w_2(~mask) = 0;
%t1w_1(~mask) = 0;
%t1w_3(~mask) = 0;

% Scale the values of t1w_2 and t1w_4
 cut = M0<(M-2*S);
%t1w_2 = t1w_2.*10000;
% t1w_2(cut) = 100;
t1w_4 = t1w_4.*10000;
M=mean(t1w_4(brainMask));
S=std(t1w_4(brainMask));

up=min(10000,M+3*S);
dwon=max(0,M-3*S);
 t1w_4(t1w_4<dwon)=dwon;
 t1w_4(t1w_4>up)=up;    

%% IV. Save out the resulting nifti files

% Set the file names for writing to disk
%Headfilefile = fullfile(outDir,'headMask.nii.gz');
%T1wfs_1file  = fullfile(outDir,'T1wfs_1.nii.gz');
%T1wfs_2file  = fullfile(outDir,'T1wfs_2.nii.gz');
%T1wfs_3file  = fullfile(outDir,'T1wfs_3.nii.gz');
T1wfs_4file  = fullfile(outDir,'T1wfs_4.nii.gz');

% Write them to disk
disp('2.  saving a syntetic T1w images');
%dtiWriteNiftiWrapper(single(mask), xform,Headfilefile);
%dtiWriteNiftiWrapper(single(t1w_1), xform,T1wfs_1file);
%dtiWriteNiftiWrapper(single(t1w_2), xform, T1wfs_2file);
%dtiWriteNiftiWrapper(single(t1w_3), xform, T1wfs_3file);
dtiWriteNiftiWrapper(single(t1w_4), xform, T1wfs_4file);

return

