function [mrQ]=mrQ_arrangeData(dataDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,channels,useNiftiFlag,mrQ)
% 
% mrQ_arrangeData(dataDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,channels,useNiftiFlag);
% 
% 
% Arrange data from a series of scans gathered at the CNI. This function
% will take a dataDir and organize the data/folders within it based upon
% the series numbers (SEIR_seriesNumbers & SPGR_seriesNumbers). It will
% dump all the data to a raw diretory at the top level of dataDir and copy
% data from there into seperate SEIR*/data and SPGR*/data directories at
% the same level as raw. The file operations here assume that each of your
% SEIR and SPGR acquisitions are in serperate folders with unique series
% numbers in their names.
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
%   SEIR_seriesNumbers- A cell array containing the SEIR_epi series numbers
%                       from dataDir/raw.  (if this input exists a SEIR dir
%                       will be created with data and fitT1_GS directories
%                       inside it). The scan dicom will be extracted from
%                       the data Dir. If SEIR_epi already exists and a new
%                       one will be created.
% 
%   SPGR_seriesNumbers- A cell array with the SPGR series numbers from
%                       dataDir/raw. This input creates a SPGR directory
%                       with data inside it. The scan dicom will be
%                       extracted from the data Dir. If SPGR already exists
%                       a new one will be created.
% 
%   channels          - A 1xN array (where N is the number of series). If
%                       any of the SPGR data has multi-channel dicoms make
%                       the series number equal to the number of channels.
%                       The default is zero.
% 
%   useNiftiFlag      - A 1xN logical array (where N is the number of
%                       series). If it is set to 1 the code will use the
%                       nii.gz file and the channel variable input is not
%                       used. 
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

if (~exist('channels','var') || isempty(channels))  && exist('SPGR_seriesNumbers','var')
    % This could be a check for the raw directory - and if it does not
    % exist then we can create it. 
    channels(1:length(SPGR_seriesNumbers)) = 0;

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


% Set the path to the raw directory, in which all data will now exist.
%mrQ.name=fullfile(dataDir,'mrQ_params');


rawDir = fullfile(dataDir,'raw');

mrQ = mrQ_Set(mrQ,'raw',rawDir); % mrQ.rawDir=rawDir; % This should use mrQ_Set
eval(['!  rm ' mrQ.RawDir '/mrQ_params.mat'])
%save(mrQ.name,'mrQ');

%% Arrange SEIR data

if ( exist('SEIR_seriesNumbers','var') && ~isempty(SEIR_seriesNumbers) )
    mrQ = mrQ_Set(mrQ,'seir',SEIR_seriesNumbers); % mrQ.SEIR_seriesNumbers=SEIR_seriesNumbers;
    
    % Initialize counters
    num = 1; ex = 0;
    
    % Make a SEIR_epi Dir. If there is one already it will make another.
    while ex == 0 
        SEIRepiDir = fullfile(dataDir,['SEIR_epi_' num2str(num)]);
        
        if ( ~exist(SEIRepiDir,'dir') )
            mkdir(SEIRepiDir);

            SEIRepiDir_dat = fullfile(SEIRepiDir,'data'); 
            mkdir(SEIRepiDir_dat);
            
            SEIRepiDir_fit = fullfile(SEIRepiDir,'fitT1_GS');
            mkdir(SEIRepiDir_fit);

            for i = 1:numel(SEIR_seriesNumbers)
                
                Dat_Dir = [rawDir '/dicoms/' SEIR_seriesNumbers{i} '*'];
                Dat_Dir = dir(Dat_Dir);
                Dat_Dir = [rawDir '/dicoms/' Dat_Dir.name];
                
                % This is duplicating data - we don't want this. 
                if exist(Dat_Dir,'dir')
                    cmd = sprintf('!cp -r %s -d %s',Dat_Dir, SEIRepiDir_dat);
                    eval(cmd);
                else
                    % Assume the directory is within a .zip file
                    seirFile = fullfile(rawDir,'/dicoms/',[SEIR_seriesNumbers{i} '*.zip']);
                    seirFile = ls(seirFile);

                    cmd = sprintf('! unzip %s -d %s',seirFile(1:end-1),SEIRepiDir_dat);
                    eval(cmd)
                end
                ex=1;
            end
        end
        num=num+1;
    end
    mrQ.SEIRepiDir=SEIRepiDir;
    mrQ.SEIR_done=0;

end
save(mrQ.name,'mrQ');


%% Arrange SPGR data

if  exist('SPGR_seriesNumbers','var') && ~isempty(SPGR_seriesNumbers)  
    mrQ = mrQ_Set(mrQ,'spgr', SPGR_seriesNumbers); %mrQ.SPGR_seriesNumbers=SPGR_seriesNumbers;
    if (~exist('channels','var') || isempty(channels))
        channels(1:length(SPGR_seriesNumbers))=0;
    end
    if (~exist('useNiftiFlag','var') || isempty(useNiftiFlag))
        useNiftiFlag(1:length(SPGR_seriesNumbers))=0;
    end
    
    % Make a SPGRDir Dir. if there is one allready it will make one more
    num = 1; ex = 0;
    
    while ex == 0 
        SPGRDir=fullfile(dataDir,['SPGR_' num2str(num)]);
        
        if ~exist(SPGRDir,'dir') 
            mkdir(SPGRDir);

            SPGRDir_dat = fullfile(SPGRDir,'data'); mkdir(SPGRDir_dat);

            for i=1:numel(SPGR_seriesNumbers)

                Niifile=dir([rawDir '/' SPGR_seriesNumbers{i} '*.nii.gz']);
                Niifile=[rawDir '/'  Niifile.name];

                Dat_Dir=[rawDir '/dicoms/' SPGR_seriesNumbers{i} '*'];
                
                Dat_Dir=dir(Dat_Dir);
                if Dat_Dir.name(end-3:end)=='.zip'
                    cmd = sprintf('! unzip %s -d %s',[rawDir '/dicoms/' Dat_Dir.name],[rawDir '/dicoms/']);
                    eval(cmd);
                    
                    Dat_Dir1=[rawDir '/dicoms/' SPGR_seriesNumbers{i} '*'];
                    Dat_Dir=dir(Dat_Dir1);
                    eval( ['! mv ' rawDir '/dicoms/' Dat_Dir(2).name ' ' rawDir '/dicoms/_' Dat_Dir(2).name])
                    Dat_Dir=Dat_Dir(1);
                end
                

                Dat_Dir=[rawDir '/dicoms/' Dat_Dir.name];
                
                if exist(Dat_Dir,'dir')
                    if ( channels(i) == 0 && useNiftiFlag(i) == 0 )
                        cmd = sprintf('!cp -r %s -d %s',Dat_Dir, SPGRDir_dat);
                        eval(cmd);
                        
                    elseif ( useNiftiFlag(i) == 1 )
                        Dirname = [rawDir '/dicoms/' SPGR_seriesNumbers{i} '*'];
                        Dirname = dir(Dirname);
                        SPGRdataDir = fullfile(SPGRDir_dat, Dirname.name);
                        mkdir(SPGRdataDir);
                        cmd = sprintf('!cp -r %s -d %s',Niifile, SPGRdataDir);
                        eval(cmd);
                        files = dir(fullfile(Dat_Dir, '*MR*'));
                        if isempty(files)
                            files = dir(fullfile(Dat_Dir, '*dcm*'));
                        end
                        
                        cmd = sprintf('!cp %s %s',fullfile(Dat_Dir,files(1).name),fullfile(SPGRdataDir,files(1).name));
                        eval(cmd);
                        
                    elseif ( channels(i) ~= 0 && useNiftiFlag(i) == 0 )
                        Dirname = [rawDir '/dicoms/' SPGR_seriesNumbers{i} '*'];
                        Dirname = dir(Dirname);
                        SPGRdataDir = fullfile(SPGRDir_dat, Dirname.name);
                        mkdir(SPGRdataDir);
                        files = dir(fullfile(Dat_Dir, '*MR*'));
                        for ii = channels(i)+1:channels(i)+1:length(files) 
                            cmd=sprintf('!cp %s %s',fullfile(Dat_Dir,files(ii).name),fullfile(SPGRdataDir,files(ii).name));
                            eval(cmd)
                        end
                    end
                else
                    % Assume that the files are in a zip file.
                    seirFile = [rawDir '/' SPGR_seriesNumbers{i} '*.zip'];
                    seirFile =  ls(seirFile);
                    if ( channels(i) == 0 )
                        cmd = sprintf('! unzip %s -d %s',seirFile(1:end-1),SPGRDir_dat); %#ok<NASGU>

                    elseif ( channels(i)~=0 )
                        tmpDir = fullfile(rawDir,'tmp');
                        mkdir(tmpDir);
                        
                        cmd = sprintf('! unzip %s -d %s',seirFile(1:end-1),tmpDir);
                        eval(cmd);
                        
                        [~, Dirname] = fileparts(seirFile(1:end-5));
                        SPGRdataDir = fullfile(SPGRDir_dat, Dirname);
                        mkdir(SPGRdataDir);

                        files = dir(fullfile(tmpDir, Dirname,'*.dcm'));

                        for ii = channels(i)+1:channels(i)+1:length(files);
                            cmd = sprintf('!cp %s %s',fullfile(tmpDir,Dirname,files(ii).name),fullfile(SPGRdataDir,files(ii).name));
                            eval(cmd)
                        end
                        rmdir(tmpDir,'s');
                    end
                end
                ex = 1;
            end
        end
        num = num+1;
    end
    mrQ.SPGR=SPGRDir;
    mrQ.SPGR_T1fit=0;
    mrQ.SPGR_init=0;
    mrQ.SPGR_PDfit=0;
    mrQ.SPGR_PDBuild=0;
    

end


%% Start a log file that will track what was done. 

%  Save out the variables and the directories so we can feed this stuff
%  into the next functions in the pipeline. 
mrQ.Arange_Date=date;

save(mrQ.name,'mrQ');
return
