function [data, extra, xform, saveName] = mrQ_initSEIR(mrQ,SEIRdir,alignFlag,ReOrder)
%
% [data, extra, xform, saveName] = mrQ_initSEIR_ver2(mrQ,SEIRdir,alignFlag)
%
% This function loads all DICOM data found within 'SEIRdir' into a single
% matrix ('data'); aligns each series to the first; and rearranges the
% 'data' matrix to generate a SEIR_Dat.mat file, containing all the DICOM
% image data in the variable 'data' (4-D double) as well as the transform
% used to align the SEIR series in the variable 'xform'. Also in the .mat
% file is the variable 'extra' which contains 'tVec' (a 1xN vector which
% contains the inversion times for all N of the series) and 'T1Vec'
% [1:5000]. Each of these variables are returned as output.
%
% INPUTS:
%       mrQ         - The mrQ structure.
%
%       SEIRdir     - The directory containing the SEIR epi data,
%                     organized in separate folders and containing the raw
%                     DICOM images from each of the SEIR EPI acquisitions.
%
%       alignFlag   - Set to 1 if you want to align the SEIR slices. 
%                   - If the SEIR data has more than a few slices, it would
%                     be a good idea to try aligning them.
%
%
%
% OUTPUTS:
%       data        - A matrix with all the DICOM data
%
%       extra:      - A structure which contains two vectors:
%                         --tVec: Inversion time for each series
%                                = [2400 1200 400 50] 
%                         --T1Vec: Possible T1 values, from 1 to 5000 msec
%                                = [1:5000]
%
%       xform       - A matrix, describing the transform computed to align
%                     the SEIR slices to the first. Returned from relaxAlignAll
%
%       saveName    - path to the saved data file
%
% These Outputs will be used as Inputs to subsequent functions by Barral et
% al., such as mrQ_fitSEIR_T1.m.
%
% USAGE NOTES:
%       Note that this function works only with DICOMs. The different
%       inversion time DICOMs need to be under a directory called "data".
%       The data directory should be under the SEIR path, then the SEIRdir.
%
%       This function is just a modification, made by Aviv Mezer in June
%       2011, on the getData.m function which was written by J. Barral, M.
%       Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009.
%
% EXAMPLE USAGE:
%       SEIRdir   = '/baseDir/QI/20110622_0582/SEIR_epi_1';
%       alignFlag = 1;
%
%       [data extra xform] = mrQ_initSEIR(SEIRdir,alignFlag);
%
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
%
% (C) Stanford University, VISTA
%

%% I. Check Inputs

if notDefined('SEIRdir') || ~exist(SEIRdir,'dir')
    SEIRdir = uigetdir(pwd,'Select your base data directory.');
end

% We could default to align if there are more than X number of series.
if notDefined('alignFlag') || isempty(alignFlag)
    alignFlag = 0;
end


% Set the path to the directory containing the data folders.
loadPath = fullfile(SEIRdir, 'data');
if ~exist(loadPath,'dir'), mkdir(loadPath); end

% Set the name of the file to which we will save the results in loadPath.
saveName = fullfile(loadPath, 'SEIR_Dat');

%% II. Load raw data into a single structure, called "d".

d=mrQ_input2struct(mrQ.inputdata_seir);

if ~notDefined ('ReOrder') 
    % ReOrder needs to be a value smaller than length(d);
    % we noticed that since the images are reistered to the first image,
    % the quality of registration may vary depending on the alignment order.
    
    dtmp=d;
    d(1)=d(ReOrder);
    d(ReOrder)=dtmp(1);
end

%% III. Align series in "d" and deal with complex data

% Align & ~Complex
if (alignFlag == 1) 
    %% NOTE:
    % Originally, we used SPM 8 rigid body registration code applied by
    % RFD; for details, see relaxAlignAll. This was tested, and it worked
    % on GE data. However, this code is not working for our Siemens data;
    % we therefore use fsl implementation (for details, see
    % mrQ_fslAlignCall). It is still to be tested whether all data need to
    % be analyzed using this code.
    
    %    SPM
    % Only gets the resolution of the first 3 dimensions. The fourth is
    %likely the TR. Then align and reslice with spm:
    
            mm = d.mmPerVox; mm = mm(1:3);
            [d, xform] = relaxAlignAll(d, [], mm, false, 1);
            
%     % FSL
%       This is another option to get alignment. It probably doesn't work
%       well with different contrasts -- check before using it!

%     Dpath=fullfile(mrQ.SEIRepiDir,'data');
%     niilist=mrQ.inputdata_seir.name;
%     [d, xform]=mrQ_fslAlignCall(Dpath,d,niilist);
%     
else
  xform=  d(1).imToScanXform;
end
%% IV. Initialize data matrix
% Determine the number of rows (nRows), number of Columns (nCol), number of
% slices (nSlice), and number of series (nSeries). Initialize the data
% matrix and extra structure.
nRow    = size(d(1).imData, 1);
nCol    = size(d(1).imData, 2);
nSlice  = size(d(1).imData, 3);
nSeries = length(d);

data       = zeros(nRow,nCol,nSlice,nSeries);
extra.tVec = zeros(1,nSeries); % One series corresponds to one TI (SEIR)
if alignFlag==0
    xform=d(1).imToScanXform;
end
% Populate 'data' with image data in 'd(k)'
for k = 1:nSeries
    dataTmp = d(k).imData;
    dataTmp = double(squeeze(dataTmp));
    
    for ss = 1:nSlice
        
        data(:,:,ss,k) = dataTmp(:,:,ss);
        
    end
    
    extra.tVec(k) = d(k).inversionTime;
end


%% V. Save the data.

extra.T1Vec = 1:5000; % This can be reduced to speed up the code

% TI = extra.tVec; % Not returned or used

save(saveName,'data','extra','xform')

return
