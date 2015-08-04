function mrQ_deleteIDfile(mrQ)

funcpath=which('mrQ_arrangeOutPutDir.m'); 
[mrQpath,~]=fileparts(funcpath); 

% mrQpath =/home/shai.berman/Documents/Code/git/mrQ

sge=fullfile(mrQpath,'sge_subjects');
filesPath=[sge,'/*'];
delCommand=sprintf('rm %s', filesPath);
[status, result]=system(delCommand);

