function mrQ_deleteIDfile(mrQ)

funcpath=which('mrQ_arrangeOutPutDir.m'); 
[mrQpath,~]=fileparts(funcpath); 
sge=fullfile(mrQpath,'sge_subjects');

% mrQpath =/home/shai.berman/Documents/Code/git/mrQ
k=mrQ.sub(isstrprop(mrQ.sub, 'digit')); 
id=str2num(k(1:8));
sgefile=num2str(id); %% change 

filesPath=fullfile(sge,['_',sgefile,'.mat']);
delCommand=sprintf('rm %s', filesPath);
[status, result]=system(delCommand);

