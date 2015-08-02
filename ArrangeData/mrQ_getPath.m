function [mrQStructPath]= mrQ_getPath(IDnum)

funcpath=which('mrQ_arrangeOutPutDir.m');
[mrQpath,~]=fileparts(funcpath);

% mrQpath =/home/shai.berman/Documents/Code/git/mrQ

directory=fullfile(mrQpath,'sge_subjects'); % some directory with all of the mrQ structures. 
%%
% load all files in library
    filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
    filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
    
    numOfFiles = length(filesInDir);
    mrQStructPath=[];
    i=1;
    while(i<=numOfFiles)
        filename = filesInDir(i).name;      % Store the name of the file
        ID=str2double(filename(2:end-4));%         take only the number- remove underscore, and .mat
       if ID==IDnum
           IDfilename=fullfile(directory, filename);
           load(IDfilename);
           mrQStructPath=(sge_info);
           return
       end
    end

if isempty(mrQStructPath)
    error('the path to the given mrQ ID was not found')
end