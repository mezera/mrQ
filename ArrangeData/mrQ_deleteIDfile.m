function mrQ_deleteIDfile(mrQ)

funcpath=which('mrQ_arrangeOutPutDir.m'); 
[mrQpath,~]=fileparts(funcpath); 

% mrQpath =/home/shai.berman/Documents/Code/git/mrQ
k=mrQ.sub(isstrprop(mrQ.sub, 'digit')); 
id=str2num(k(1:8));
sgefile=num2str(id); %% change 

sge=fullfile(mrQpath,'sge_subjects');
filesPath=[sge,'_',sgefile];
delCommand=sprintf('rm %s', filesPath);
[status, result]=system(delCommand);

