function mrQ_createIDfile(mrQ)

funcpath=which('mrQ_arrangeOutPutDir.m'); 
[mrQpath,~]=fileparts(funcpath); 

% mrQpath =/home/shai.berman/Documents/Code/git/mrQ

sge=fullfile(mrQpath,'sge_subjects');
if ~exist(sge, 'dir')
    mkdir(sge) 
end

k=mrQ.sub(isstrprop(mrQ.sub, 'digit')); 
id=str2num(k(1:8));
sgefile=num2str(id); %% change 

mysgefile=fullfile(sge,strcat('_',sgefile));
sge_info=mrQ.name; %location of the mrQ structure (.mat)
save(mysgefile, 'sge_info');
