function [data, extra, xform, saveName] = mrQ_initSEIR_ver2(mrQ,SEIRdir,alignFlag)
%
% [data extra xform saveName] = mrQ_initSEIR(mrQ,SEIRdir,alignFlag)
%
% Loads all dicom data found within 'SEIRdir' into a single matrix
% ('data'), aligns each series to the first, rearranges the 'data' matrix
% and generates a SEIR_Dat.mat file containing all the DICOM image data in
% the variable 'data' (4-D double) as well as the transform used to align
% the SEIR series in the variable 'xform'. Also in the .mat file is the
% variable 'extra' which contains 'tVec' - (1xN) the inversion times for
% each series (N) and 'T1Vec' - [1:5000]. Each of these variables are
% returned.
%
% INPUTS:
%       SEIRdir     - The directory containing your SEIR epi data,
%                     organized in seperate folders containg the raw dicom
%                     images from each of the SEIR EPI acquisitions.
%
%       alignFlag   - Set to 1 if you want to align the SEIR slices. - If
%                     the SEIR data has more than a few slices it is a good
%                     idea to try to align them.
%

%
% OUTPUTS:
%                     ** (Check with Aviv to make sure this is correct) **
%       data        - matrix with all the dicom data
%
%       extra:      - tVec:  [2400 1200 400 50] = Inversion time for each
%                                                 series.
%                   - T1Vec: [1x5000 double] = ?
%
%       xform       - The transform computed to align the SEIR slices to
%                     the first. Returned from relaxAlignAll
%
%       saveName    - path to the saved data file
%
% USAGE NOTES:
%       Note that this function works only with dicoms. The different
%       inversion time dicoms need to be under a directory called "data".
%       The data directory should be under the SEIR path, then the SEIRdir.
%       this function is just a  modification made by aviv mazer in june
%       2011 on  the getData.m that was written by J. Barral, M.
%       Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%
% EXAMPLE USAGE:
%       SEIRdir   = '/baseDir/QI/20110622_0582/SEIR_epi_1';
%       alignFlag = 1;
%
%       [data extra xform] = mrQ_initSEIR(SEIRdir,alignFlag);
%
%
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
%
% (C) Stanford University, VISTA
%

%% Check INPUTS

if notDefined('SEIRdir') || ~exist(SEIRdir,'dir')
    SEIRdir = uigetdir(pwd,'Select your base data directory');
end

% We could default to align if there are more than X number of series.
if notDefined('alignFlag') || isempty(alignFlag)
    alignFlag = 0;
end



% Set the path to the directory containing the data folders
loadPath = fullfile(SEIRdir, 'data');
if ~exist(loadPath,'dir'), mkdir(loadPath); end

% Set the name of the file that we will save the results to in loadPath.
saveName = fullfile(loadPath, 'SEIR_Dat');


%% Load raw  data into a single structure (d)

d=mrQ_input2Stuck(mrQ.inputdata_seir);


%% Align series in d and deal with complex data

% Align & ~Complex
if (alignFlag == 1) 
    %% NOTES
    % originatly we used spm 8 ridge body registration code applied
    % by RFD for detail see relaxAlignAll. this was tested and work on
    % GE data.
    % this code is not working for our siemens data we there for use
    % fsl implemntaion for details see mrQ_fslAlignCall. it still to test if all data need to be
    % use this code.
    %
    
    %    SPM
    %Only get the resolution of the first 3 dimensions. The fourth is
    %likely the TR. The align and reslice with spm
            mm = d.mmPerVox; mm = mm(1:3);
            [d xform] = relaxAlignAll(d, [], mm, false, 1);
            
%     %fsl
%       this another option to get alignment. it probably doesn't work
%       well with different contrasts. check before using it!

%     Dpath=fullfile(mrQ.SEIRepiDir,'data');
%     niilist=mrQ.inputdata_seir.name;
%     [d, xform]=mrQ_fslAlignCall(Dpath,d,niilist);
%     
    
end
% Determine number of rows (nRows) number of Columns (nCol) and number of
% slices (nSlice) and number of series (nSeries) and initialize the data
% matrix and extra structure.
nRow    = size(d(1).imData, 1);
nCol    = size(d(1).imData, 2);
nSlice  = size(d(1).imData, 3);
nSeries = length(d);

data       = zeros(nRow,nCol,nSlice,nSeries);
extra.tVec = zeros(1,nSeries); % One series corresponds to one TI (SEIR)

% Populate 'data' with image data in 'd(k)'
for k = 1:nSeries
    dataTmp = d(k).imData;
    dataTmp = double(squeeze(dataTmp));
    
    for ss = 1:nSlice
        
        data(:,:,ss,k) = dataTmp(:,:,ss);
        
    end
    
    extra.tVec(k) = d(k).inversionTime;
end


%% Save the data out to a file

extra.T1Vec = 1:5000; % This can be reduced to speed up the code

% TI = extra.tVec; % Not returned or used

save(saveName,'data','extra','xform')


return
