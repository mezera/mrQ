function [curentTime]=mrQ_cheachSGE(curentTime,timeLimit,opt,SGEEval,sgwoutPath)
% %
% %
% %
% %
% %
% 
% % If too much time has elapsed then we recall the grid;
%          if curentTime >timeLimit  % 24hours
%              ch=[1:jumpindex:length(opt.wh)]; %the nude filre name
%              k=0;
%              reval=[]
%              for ii=1:length(ch),
% %                 
%                  ex=['_' num2str(ch(ii)) '_'];
%                 if length(regexp(list, ex))==0,
%                      k=k+1;
%                      reval(k)=(ii); % we make a list of the grid run that are not done yet
% %                     
%                  end
%              end;
%              
%              
%              for i=1:length(reval)
%                              [status, result] =system([' find ' location '-type f -mmin -300']);
%              [status, result]=system(' find ~/sgeoutput/*   '-type f -mmin -15')
%              
%              
%              end
%              eval(opt.clean);
%              evan(opt.gridRun)
%              t=0;
%              