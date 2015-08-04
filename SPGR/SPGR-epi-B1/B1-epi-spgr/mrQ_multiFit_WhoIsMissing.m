function MissingFileNumber=mrQ_multiFit_WhoIsMissing( dirname,N_Vox2Fit,jumpindex)
% MissingFileNumber=mrQ_multiFit_HowIsMissing( dirname,N_Vox2Fit,jumpindex)
%

%%
MissingFileNumber = [];
%make a list of files
list  = ls(dirname);

%% make a list of expected files
ch    = 1:jumpindex:N_Vox2Fit;

k     = 0; % counting integer
%% check if all files are there
for ii=1:length(ch), % run over expected files
    ex=['_' num2str(ch(ii)) '_']; % a string that the saved file needs to have
    if length(regexp(list, ex))==0, % check if the such file exists
        k=k+1;
        MissingFileNumber(k)=(ii); % identify the files that are missing
    end
end
