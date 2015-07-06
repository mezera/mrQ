function GridFit_done=mrQ_Gridcheck(opt_Log_name,SunGrid,CallType)
%%
% This loop checks if all the outputs have been saved or waits until
% they are all done, if needed it will run again job that are missing

if (~exist(opt_Log_name,'file')),
    
    disp(['cant find file : ' opt_Log_name  ])
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
        % check if ther are jobs on the sun grid que list
        qStatCommand    = [' qstat | grep -i  job_' sgename(1:6)];
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
