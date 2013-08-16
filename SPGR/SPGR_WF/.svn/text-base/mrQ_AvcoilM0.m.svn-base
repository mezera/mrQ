function saveFileName = mrQ_AvcoilM0(niifile,datDir)
% 
% saveFileName = mrQ_AvcoilM0(niifile,datDir)
% 
% Combine aligned MO data from each coil with the multi channel data and
% save that data out in 'fullfile(datDir, 'AligncombineCoilsM0')'.
% 
% Adds the data fields of each of the nifti files (of number N) and
% divides each element in the array by N. Saves out as nifti.
% 
% INPUTS:
%   niifile     - a cell array of nifti files that contain data from each
%                 channel (N nifti files for N channels).
%   datDir      - The directory containing the nifti files - also the
%                 directory containing 'dat_aligned.mat' - the aligned
%                 data
% 
% OUTPUT:       
%   saveFileName - The full path to the nifti file containing the combined
%                  data.
%         
% 
% SEE ALSO:
%   mrQ_multicoilM0.m
% 
% 
% 
% (C) Stanford University, VISTA 
% 
% 

% Note: This funciton should be renamed to reflect it's actual function.
%       mrQ_combineChannelData
% 


%% Combine channel data from nifti files

for i=1:length(niifile)
    P1 = readFileNifti(niifile{i});
    if i==1;
        xform = P1.qto_xyz;
        combine = zeros(size(P1.data));
    end
    combine = combine + P1.data;
end

combine = combine./i;

saveFileName = fullfile(datDir, 'AligncombineCoilsM0');

dtiWriteNiftiWrapper(single(combine), xform, saveFileName);

return