function combineFile = mrQ_multicoilM0(datDir,T1file,B1file,niifile,flip_Angles,mrQ)
% 
% combineFile = mrQ_multicoilM0(datDir,T1file,B1file,niifile,flip_Angles,sics)
% 
% Align individual channels to the reference image and perform the fit for
% each channel/coil. Calculate M0 and combine the data from each channel
% and save to a file (combineFile). #Do minimization for PD (see
% fitPD_andGain folder).
% 
% INPUTS:
%   datDir      - Directory containing 'dat_aligned.mat' - the aligned
%                 data.
%   T1file      - Full path to 'T1_lsq_GLr.nii.gz'
%   B1file      - Full path to 'B1_fit_lregGx3.nii.gz'
%   niifile     - A cellstr array of nifti file names for each of the F_A 
%                 acquisitions 
%   flip_Angles - An array of flip angles in the order that the data should
%                 be combined.

% 
% 
% OUTPUTS:
%   combineFile - File name of the nifti file containing the combined
%                 channel data
% 
% SEE ALSO:
%   mrQ_AvcoilM0.m , relaxAlignAll_multichanels.m
% 
% 
% (C) Stanford University, VISTA Lab
% 
% 


%% Check and Load INPUTS 

outFile = fullfile(datDir,'dat_aligned.mat');
disp(['Loading aligned data from ' outFile '...']);
load(outFile);
 

% Read the flip angles and TR from the actual data
flipAngles = [s(:).flipAngle];
TR         = [s(:).TR];
for ii=1:length(flipAngles)
SZraw(:,:,:,ii)=size(s(ii).imData);
end

% Load the T1 and scale the T1 values
if (~exist('T1file','var') || isempty(T1file)),
    T1file= fullfile(datDir,'T1_map_lsq.nii.gz');
end

if ~exist(T1file,'file') && isfield(mrQ,'outDir') 
    T1file=mrQ_getT1file(mrQ);
    %T1file = fullfile(mrQ.outDir,'T1_map_lsq.nii.gz');
end
T1 = readFileNifti(T1file);
T1 = double(T1.data);T1=T1.*1000;


% Load the B1
if (~exist('B1file','var') || isempty(B1file)),
    B1file= fullfile(datDir,'B1_Map.nii.gz');
end

if ~exist(B1file,'file') && isfield(mrQ,'outDir') 
    B1file = fullfile(mrQ.outDir,'B1_Map.nii.gz');
end

B1 = readFileNifti(B1file);
B1 = double(B1.data);


%% Align the channels to the reference image and combine them


% For each flip angle align the individual channel data to the reference
% image and combine them
imN=0;
FAU= unique(flipAngles);

for j=1:numel( FAU)
    
    kkk=0;
    % find all images with flip anlge j
    kk = find(flipAngles == FAU(j)); % won't kk always = j? *** originaly it was not the order of the flipangle as an input does not have to be the order of flipangle images saved in the s structure.
  %  ref   = fullfile(datDir,['Align' num2str(flipAngles(j)) 'deg']);
    
    % Loop over images with the same flip angle
    for fa = 1:length(kk)
         imN=imN+1;

        kkk=kkk+1; % Count
        nk=kk(kkk);
    
        ref   = fullfile(datDir,['Align' num2str(FAU(j)) 'deg_' num2str(kkk)]);
    mrQ_makeNiftiFromStruct(s(nk),ref,xform);
    refIM = readFileNifti(ref);
    
    % niifile here is the raw image???
    s1 = makeStructFromNifti(niifile{nk},-1,[],mrQ.permution);

    channels = length(s1)-1;
    [s11,xform1 s1] = relaxAlignAll_multichanels(s1(channels+1), refIM, mmPerVox, true, 1,s1(1:channels));
    
    
    szref = SZraw(:,:,:,nk)';
    szdat = size(s11(1).imData);
    clear refIM s11
    % Dimension checks 
    if szref(1)>szdat(1);   for k=1:length(s1); s1(k).imData(szdat(1):szref(1),:,:)=0;end;   end
    if szref(2)>szdat(2);   for k=1:length(s1); s1(k).imData(:,szdat(2):szref(2),:)=0;end;   end
    if szref(3)>szdat(3);   for k=1:length(s1); s1(k).imData(:,:,szdat(3):szref(3))=0;end;   end
    
    
    fa  = flipAngles(nk).*B1;
    fa  = fa./180.*pi;
    tr  = TR(nk);
    
    M0c = zeros(szref(1),szref(2),szref(3),channels);
    
    
    % Calculate M0 for each channel
    for i = 1:channels
        M0c(:,:,:,i) = s1(i).imData(:,:,:) ./ ( (1-exp(-tr./T1)).*sin(fa) ./ (1-exp(-tr./T1).*cos(fa)) );
    end
    
    clear s1 s11 refIM fa szref szdat
    
    % Save out the data for each channel
    SaveFilename = [ref 'MC_m0'];
    fprintf('save file: %s',SaveFilename);
    dtiWriteNiftiWrapper(single(M0c),xform, SaveFilename);
    
    clear M0c ref 
    
    Files{nk} = SaveFilename;
end
end
% Combine the data for each coil using the multi channel data and save it to
% disk
 combineFile = mrQ_AvcoilM0(Files,datDir);
 
 return




    
