function [mrQ]=mrQ_arrangeDataMT(dataDir,arrangeRawFlag,MT_seriesNumbers)
% 
% mrQ_arrangeData(dataDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,channels,useNiftiFlag);
% 
% 
% Arrange MT data from a series of scans gathered at the CNI. This function
% will take a dataDir and organize the data/folders within it based upon
% the series numbers (MT_seriesNumbers ). It will
% dump all the data to a raw diretory at the top level of dataDir and copy
% data from there into seperate MT*/data directorie at
% the same level as raw. The file operations here assume that each of your
% MT acquisitions are in serperate folders with unique series
% numbers in their names and there is a Niftii and dicom directory of the data.
% 
% 
% INPUT VARIABLES:
%   dataDir           - The path of the directory containing the raw MR
%                       data. the function assume the data are arange in
%                       the way the Stanford CNI. the data should be saved
%                       as niifti files and a dicom directory with dicom
%                       are also exsist. 
% 
%   arangeRawFlag     - Needs to be 1 to move every file from the original
%                       magnet directory (dataDir) to a raw directory
%                       inside it (if this was done set = 0, or empty).
% 
%   MT_seriesNumbers- A cell array containing the MT series numbers
%                       from dataDir/raw.  (if this input exists a MT dir
%                       will be created with data 
% 

%
% OUTPOUTs   
%                        mrQ information structure will be saved. the raw data will be move to raw Dir.
%                       the SPGR will be saved in SPGR dir.
%                       the SEIR in SEIR DIR. 
%                       each data set will have it's own Di structure for the anaysis.
  %                       
% 
% EXAMPLE USAGE:
%     dataDir = '/biac2/wandell2/data/WMDevo/adult/QuantitativeImaging/109_AL/20110622_0582';
%     arrangeRawFlag = 1;
%     SEIR_seriesNumbers = {'05' '06' '07' '08'};
%     SPGR_seriesNumbers = {'09' '10' '11' '12'};
%     channels = [32 32 32 32];
%     useNiftiFlag = [1 1 1 1]; 
% 
%     mrQ_arrangeData(dataDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,channels,useNiftiFlag);
% 
% 
% (C) Stanford, VISTA [2011]
% 
% 

%% Check INPUT variables

if notDefined('dataDir') || ~exist(dataDir,'dir')
    dataDir = uigetdir(pwd,'Select your base data directory');
end

if ~exist('arrangeRawFlag','var') || isempty(arrangeRawFlag) 
    % This could be a check for the raw directory - and if it does not
    % exist then we can create it. 
    arrangeRawFlag = 0;
end

% Make a raw directory and move all data into this folder, which
% will then be placed in dataDir. 
if arrangeRawFlag 
    oneup  = fileparts(dataDir);
    rawDir = fullfile(oneup,'raw');
    if ~isdir(fullfile(dataDir,'raw')) 
        mkdir(rawDir); 
        movefile([dataDir '/*'], rawDir);
        movefile(rawDir, dataDir);
    end
end
rawDir = fullfile(dataDir,'raw');


if ( exist('MT_seriesNumbers','var') && ~isempty(MT_seriesNumbers) )
    mrQ.MT_seriesNumbers=MT_seriesNumbers;

% Initialize counters
    num = 1; ex = 0;
    
    % Make a MT Dir. If there is one already it will make another.
    while ex == 0 
        MTDir = fullfile(dataDir,['MT_' num2str(num)]);
        
        if ( ~exist(MTDir,'dir') )
            mkdir(MTDir);
        
            MTDir_dat = fullfile(MTDir,'data'); 
            mkdir(MTDir_dat);
            
          

%% Arrange MT data


            for i = 1:numel(MT_seriesNumbers)
                Niifile=dir([rawDir '/' MT_seriesNumbers{i} '*.nii.gz']);
                Niifile=[rawDir '/'  Niifile.name];

                Dat_Dir=[rawDir '/dicoms/' MT_seriesNumbers{i} '*'];
                Dat_Dir=dir(Dat_Dir);
                Dat_Dir=[rawDir '/dicoms/' Dat_Dir.name];
                
               
                 Dirname = [rawDir '/dicoms/' MT_seriesNumbers{i} '*'];
                        Dirname = dir(Dirname);
                       MTdataDir = fullfile(MTDir_dat, Dirname.name);
                        mkdir(MTdataDir);
                        cmd = sprintf('!cp -r %s -d %s',Niifile, MTdataDir);
                        eval(cmd);
                        files = dir(fullfile(Dat_Dir, '*MR*'));
                        cmd = sprintf('!cp %s %s',fullfile(Dat_Dir,files(1).name),fullfile(MTdataDir,files(1).name));
                        eval(cmd);
              
            end
            ex=1;
        end
        num = num+1;

    end
    
    mrQ.MTDir=MTDir;
    mrQ.MT_T1fit=0;
    mrQ.MT_init=0;
    
    

end
                
%% Start a log file that will track what was done. 

%  Save out the variables and the directories so we can feed this stuff
%  into the next functions in the pipeline. 
mrQ.MT_Arange_Date=date;


return
