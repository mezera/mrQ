function GridFit_done=mrQ_Gridcheck(opt_Log_name,SunGrid,CallType,GridOutputDir)
% function GridFit_done=mrQ_Gridcheck(opt_Log_name,SunGrid,CallType)
%
% This function employs a "while" loop to check whether all the outputs
% have been saved, and will wait until they are all done. If needed, the
% function will re-run for jobs that are missing.
%
%
%  ~INPUTS~
%       opt_Log_name:    Location of the "opt" structure
%            SunGrid:    Whether to use SunGrid (default is 0, "no")
%           CallType:    Determines what type of function will be
%                           performed, assuming RunSelectedJobs is true. It
%                           will be performed in the SunGrid. Enter 1 for
%                           mrQ_fitM0boxesCall_Multi.m (default), 2 for
%                           mrQ_fitB1boxesCall.m (currently
%                           non-operational), or 3 for mrQ_fitB1LR_Call.m.
%      GridOutputDir:    The 
%
%  ~OUTPUTS~
%          GridFit_done: A true/false logical, which indicates whether the
%                              function has run its course 
%
%(C) Mezer lab, the Hebrew University of Jerusalem, Israel
%  2015

%% I. Load opt file
if (~exist(opt_Log_name,'file')),
    
    disp(['Cannot find file : ' opt_Log_name  ])
    error
else
    load  (opt_Log_name)
end


if notDefined('SunGrid');SunGrid=0;end
if notDefined('CallType');CallType=1;end
if notDefined('GridOutputDir')
    GridOutputDir=pwd;
end
GridFit_done=false;
%%
if CallType==1 || CallType==2
    fNum=ceil(length(opt.wh)/opt.jumpindex);
elseif CallType==3
    fNum=ceil(opt.N_Vox2Fit/opt.jumpindex);
end

sgename=opt.SGE;
fullID=sgename(isstrprop(sgename, 'digit'));
jobname=(fullID(1:3));

tic
while GridFit_done~=true
    
    % List all the files that have been created from the call to the SGE
    list=ls(opt.dirname);
    
    % Check if all the files have been made.  If they are, then collect
    % all the nodes and move on.
    if length(regexp(list, '.mat'))>=fNum,
        GridFit_done=true;
        
        % Once we have collected all the nodes, we delete the SGE output
        eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        
    else
        % Check if there are jobs in the SGE queue
        qStatCommand    = [' qstat | grep -i  job_' jobname];
        [status result] = system(qStatCommand);
        tt=toc;
        if (isempty(result) && tt>60)
            % Check if 1 minute has passed. If there are no jobs waiting to
            % be finished (and we don't have all the jobdone) then we will
            % need to re-run it.
            
            RunSelectedJob=true;
            if CallType==1
                if isfield(opt,'Reg') % Allow different fit methods
                    mrQ_fitM0boxesCall_Multi(opt_Log_name,SunGrid,RunSelectedJob)
                else
                    mrQ_fitM0boxesCall(opt_Log_name,SunGrid,RunSelectedJob)
                end
            elseif CallType==2
                mrQ_fitB1boxesCall(opt_Log_name,SunGrid,RunSelectedJob);
            elseif CallType==3
                mrQ_fitB1LR_Call(opt_Log_name,SunGrid,RunSelectedJob);
            end
        end
        
    end
end

% If the job was finished, remove all SGE outputs. 
if  GridFit_done
    jobname=(fullID(1:3));
    filesPath=[GridOutputDir,'/job_',num2str(jobname),'*'];
    delCommand=sprintf('rm %s', filesPath);
    [status, result]=system(delCommand);
end

