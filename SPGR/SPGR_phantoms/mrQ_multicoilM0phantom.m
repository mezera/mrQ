function combineFile = mrQ_multicoilM0phantom(datDir,T1file,B1file,niifile,flip_Angles,mrQ)
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

outFile = fullfile(datDir,'dat_aligned.mat');
disp(['Loading aligned data from ' outFile '...']);
load(outFile);
 

% Read the flip angles and TR from the actual data
flipAngles = [s(:).flipAngle];
TR         = [s(:).TR];
clear s;




% Load the T1 and scale the T1 values
if (~exist('T1file','var') || isempty(T1file)),
    T1file= fullfile(datDir,'T1_map_lsq.nii.gz');
end
T1 = readFileNifti(T1file);
T1 = double(T1.data);T1=T1.*1000;


% Load the B1
if (~exist('B1file','var') || isempty(B1file)),
    B1file= fullfile(datDir,'B1_Map.nii.gz');
end
B1 = readFileNifti(B1file);
B1 = double(B1.data);


%% Align the channels to the reference image and combine them
kkk=1;

% For each flip angle align the individual channel data to the reference
% image and combine them
for j=1:length(flip_Angles)
    kk = find(flipAngles==flip_Angles(j));
    
     
    % The ref image here will be an indivudual nifti file containing the
    % data from a given flip angle. 
    ref = fullfile(datDir,['Align' num2str(flip_Angles(j)) 'deg']);
    
    if length(kk)>1
        
        kk=kk(kkk);
       ref   = fullfile(datDir,['Align' num2str(flipAngles(j)) 'deg_' num2str(kkk)]); 
        kkk=kkk+1;
    end
  
    % niifile here is the raw image???
    s1 = makeStructFromNifti(niifile{j},-1,[],mrQ.permution);
    
   channels = length(s1)-1;
    
    szref = size(s1(1).imData);
    
    
    fa  = flipAngles(kk).*B1;
    fa  = fa./180.*pi;
    tr  = TR(kk);
    
    M0c = zeros(szref(1),szref(2),szref(3),channels);
    
    
    % Calculate M0 for each channel
    for i = 1:channels
        M0c(:,:,:,i) = double( s1(i).imData(:,:,:)) ./ ( (1-exp(-tr./T1)).*sin(fa) ./ (1-exp(-tr./T1).*cos(fa)) );
    end
    
    clear s1 s11 refIM fa szref szdat
    
    % Save out the data for each channel
    SaveFilename = [ref 'MC_m0'];
    fprintf('save file: %s',SaveFilename);
    dtiWriteNiftiWrapper(single(M0c),xform, SaveFilename);
    
    clear M0c ref kk
    
    Files{j} = SaveFilename;
end

% Combine the data for each coil using the multi channel data and save it to
% disk
 combineFile = mrQ_AvcoilM0(Files,datDir);
 
 return




    
