function [s2,mmPerVox,xform] = mrQ_multicoilWeighting_phantoms(datDir,niiFile,flipAngles)
% 
% [s2,mmPerVox,xform] = mrQ_multicoilWeighting(datDir,niiFile,flipAngles)
% 
% Determine coil weighting given multicoil SPGR data.
% code). 
% 
% INPUTS:
%       datDir -    Contains dat_aligned.mat - which is the result of
%                   running mrQ_initSPGR.m (this function is now called
%                   from within that
% 
%       niiFile -   A cell array containing the path to each of the nifti
%                   files produced in mrQ_initSPGR.
% 
%       flipAngles- A 1xN vector with one entry for each of the flip
%                   angles. This is also in dat_aligned.mat (in
%                   s(:).flipAngle).
% OUTPUTS:
%       s2       - xformed/aligned data structure containing only data from
%                  the top 8 channels (8 set by opt.numIn). Returned from
%                  mrQCalweights. * Could use more explanation
% 
%       mmPerVox - Millimeter per voxel in 'inFile'.
% 
%       xform    - The transform in dat_aligned.mat ('inFile')
% 
% 
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
% 
% 
% (C) Stanford University, VISTA Lab
% 
% 

%% Check INPUTS and load data

if notDefined('datDir') || ~exist(datDir,'dir')
    datDir = uigetdir(pwd,'Select your aligned SPGR data directory');
end

if notDefined('niiFile')
    error('You must provide a cell array of nifti files to be processed');
end

% Load the data from the mat file.
inFile = fullfile(datDir,'dat_aligned.mat');
disp(['Loading aligned data from ' inFile '...']);
load(inFile);

% Get flip angles for each acquisition
if notDefined('flipAngles'),
    flipAngles = [s(:).flipAngle];
end

% Repitition times for each acquisition *** Not sure we need this.
TR = [s(:).TR];


%% Set up the grids that will be used to determine which channel is best *** 

% Dimension handling
sz   = size(s(1).imData);
boxS = round(24./mmPerVox);
even = find(mod(boxS,2)==0);

% Make sure that boxS does not have any even # dimension
boxS(even) = boxS(even)+1;
opt.boxS   = boxS;

% Determine the center of the boxes we will fit and create the grids and
% set grid options for the coil weighting
[opt.X,opt.Y,opt.Z] = meshgrid(round(boxS(1)./2):boxS(1):sz(1) , ... 
    round(boxS(2)./2):boxS(2):sz(2) , round(boxS(3)./2):boxS(3):sz(3));
opt.donemask = zeros(size(opt.X));
opt.HboxS    = (boxS-1)/2;
opt.numIn    = 8; % The number of top performing coils to keep ***
opt.sz       = sz;


%% Align data and calculate weights

for j=1:numel(flipAngles)
    
    kk = find(flipAngles == flipAngles(j)); % won't kk always = j? ***

    % Create a reference nifti file (refIM) from the data structure (s)
    ref   = fullfile(datDir,['Align' num2str(flipAngles(j)) 'deg']);
    mrQ_makeNiftiFromStruct(s(kk),ref,xform); 
    refIM = readFileNifti(ref);

    
    % Get the data from the nifti into a struct- (s1) Now the data will be
    % in a structure with separate entries for each channel (-1 flag).
    s1 = makeStructFromNifti(niiFile{j},-1,[]);
    
    
    % The number of channels in the coil (It comes up +1 why?)
    channels = length(s1)-1; 
    
    
   
    
    % Size of the data for a given flip angle volume
    szref = size(s(kk).imData);
    % Size of the data structure returned form relaxAlignAll_multichannels
    szdat = size(s1(1).imData);
    
    
    % *** The data structures are different sizes. Handle that here.
    if szref(1)>szdat(1); for k=1:length(s1); s1(k).imData(szdat(1):szref(1),:,:)=0; end; end
    if szref(2)>szdat(2); for k=1:length(s1); s1(k).imData(:,szdat(2):szref(2),:)=0; end; end
    if szref(3)>szdat(3); for k=1:length(s1); s1(k).imData(:,:,szdat(3):szref(3))=0; end; end
    
    
    % Do the coil weighting and return the weighted data in s2. For the
    % first flip angle (j==1). For the rest of the volumes (j>1) Incoil
    % will be passed through. Incoil is a 8x1 vector of indices of sorted
    % channel means in descending order of the 8 coils (opt.numIn) with the
    % highest mean signal assigined within mrQ_calculateCoilWeights
    if j==1,
        [s2(j) Incoil] = mrQ_calculateCoilWeights(s1,opt,s(kk),[]);
    else
        [s2(j) dd]     = mrQ_calculateCoilWeights(s1,opt,s(kk),Incoil);
    end

end


%% Save out the data based on coil weights

outFile = fullfile(datDir,'dat_alignedBest.mat');
disp(['Saving aligned data to: ' outFile]);

% mmPerVox and xform are loaded into the workspace from 'inFile'
save(outFile,'s2','mmPerVox','xform'); % *** Should this save out xform1? 

return
