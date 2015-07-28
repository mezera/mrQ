function [mrQpath]= mrQ_getPath(IDnum)

directory=fullfile(pwd,'sge'); % some directory with all of the mrQ structures. 
load(fullfile(directory,'sge_info'));
numOfFiles = length(sge_info);
mrQpath=[];

for i=1:numOfFiles
    if sge_info(i).ID==IDnum
        mrQpath=sge_info(i).path;
        return
    end
end

if isempty(mrQpath)
    error('the path to the given mrQ ID was not found')
end