function [s,xform,mmPerVox,niiFiles,flipAngles,mrQ] =  mrQ_initSPGR_pantom(spgrDir,mrQ)
%

if notDefined('spgrDir') || ~exist(spgrDir,'dir')
    spgrDir = uigetdir(pwd,'Select your SPGR directory');
end

dicomDir = fullfile(spgrDir,'data');



% Get the path to each of the dicom directories (d)
d = genpath(dicomDir);
if(isempty(d)), error(['Dicom dir "' d '" not found or empty.']); end
if(isunix), d = explode(':',genpath(dicomDir));
else d = explode(';',genpath(dicomDir)); end


% The first and last entries in d are the root and an empty matrix
% (resp), so we get rid of those so that d is now simply the
% directories contating the dicom/niftis.
d = d(2:end-1);

% Loop over 'd' to get the nifi paths DO THIS EARLIER - OUTSIDE OF THIS
% STATEMENT
for ii = 1:numel(d)
    ni = dir(fullfile(d{ii},'*nii.gz'));
    niiFiles{ii} = fullfile(d{ii},ni.name);
    % 'niiFiles' now contains all the paths to the nifti files and can
    % be passed into mrQ_multicoil_Weights
end





%% Load data from niftis or dicoms - reshape and permute data (nifti)



% Make a dummy structure so we will have the dicom info
s = dicomLoadAllSeries(dicomDir);

% Loop over niiFiles to get the data from the nifti and combine
% with the dicom info - reshape and permute.
for ii = 1:numel(niiFiles)
    [ s1(ii)  mrQ.coilNum(ii)]= makeStructFromNifti(niiFiles{ii},-2,s(ii),mrQ.permutation);
    s1(ii).imData=double(s1(ii).imData);
    
    % consider add those for canonical corditate space
   % canXform = mrAnatComputeCannonicalXformFromDicomXform(s1(ii).imToScanXform,size(s1(ii).data));
   % [s1(ii).imData,s1(ii).mmPerVox] = applyCannonicalXform(s1(ii).data, canXform, s1(ii).mmPerVox, false);
   % s1(ii).imToScanXform= inv(canXform*inv(s1(ii).imToScanXform));

    
end

clear s
s = s1;
clear s1


 

mmPerVox = s(1).mmPerVox(1:3);
xform=s(1).imToScanXform;

%      for i=1:numel(s)
%       bb = mrAnatXformCoords(xform,[1 1 1;size(s(i).imData)]);
%     s(i).imData = mrAnatResliceSpm(double((s(i).imData),inv(xform),bb,mmPerVox,1);
%      end

outDir = fullfile(spgrDir,['Align_'  num2str(mmPerVox(1)) '_' num2str(mmPerVox(2)) '_'  num2str(mmPerVox(3))]);
if(~exist(outDir,'dir')), mkdir(outDir); end
outFile = fullfile(outDir,'dat_aligned.mat');
save(outFile,'s', 'xform', 'mmPerVox');

for f = 1:numel(s), flipAngles(f) = s(f).flipAngle; end
fprintf('\nDetermining optimal coil weighting...\n');
% Should this return the new structure with the weighting applied?

%% i need to this in manual or no alignment way ...
%     mrQ_multicoilWeighting_phantoms(outDir,niiFiles,flipAngles,permutation);


mrQ.spgr_initDir=outDir;
mrQ.SPGR_niiFile=niiFiles;
mrQ.SPGR_niiFile_FA=flipAngles;
mrQ.SPGR_init_mmPerVox=mmPerVox;

fprintf('mrQ_initSPGR_ver2.m - COMPLETE!');

return





