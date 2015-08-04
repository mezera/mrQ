function GridFit_done=mrQ_Gridcheck(opt_Log_name,SunGrid,CallType)
% function GridFit_done=mrQ_Gridcheck(opt_Log_name,SunGrid,CallType)
%
% This function employs a "while" loop to check whether all the outputs
% have been saved, and will wait until they are all done. If needed, the
% function will re-run for jobs that are missing.
%
%
%  ~INPUTS~
%          opt_Log_name: Location of the opt file
%               SunGrid: Whether to use SunGrid (default is 0, "no")
%              CallType: Determines what type of function will be
%                        performed, assuming RunSelectedJobs is true. It
%                        will be performed in the SunGrid. Enter 1 for
%                        mrQ_fitM0boxesCall_Multi.m (default), 2 for
%                        mrQ_fitB1boxesCall.m (currently non-operational), 
%                        or 3 for mrQ_fitB1LR_Call.m.
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

GridFit_done=false;
%%
if CallType==1 || CallType==2
    fNum=ceil(length(opt.wh)/opt.jumpindex);
elseif CallType==3
    fNum=ceil(opt.N_Vox2Fit/opt.jumpindex);
end

sgename=opt.SGE;
fullID=sgename(isstrprop(sgename, 'digit'));

tic
while GridFit_done~=true
    % List all the files that have been created from the call to the
    % grid
    list=ls(opt.dirname);
    % Check if all the files have been made.  If they are, then collect
    % all the nodes and move on.
    if length(regexp(list, '.mat'))>=fNum,
        GridFit_done=true;
        % Once we have collected all the nodes we delete the sge output
        eval(['!rm -f ~/sgeoutput/*' sgename '*'])
    else
        % check if there are jobs on the sun grid queue list
        jobname=(fullID(1:3));
        qStatCommand    = [' qstat | grep -i  job_' jobname];
        [status result] = system(qStatCommand);
        tt=toc;
        if (isempty(result) && tt>60)
            % check if 1 min pass and there are no job waiting to be finish (and we don't have all the jobdone)
            %then we will need to re run it.
            
            RunSelectedJob=true;
            if CallType==1
                if isfield(opt,'Reg') % ALWO DIFFERNT FITS METHODS
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
