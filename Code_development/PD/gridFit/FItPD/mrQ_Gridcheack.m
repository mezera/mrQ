function PDGridFit_done=mrQ_Gridcheack(opt_Log_name,SunGrid,proclus)
%%
% This loop checks if all the outputs have been saved or waits until
% they are all done, it it's too long run again

if (~exist(opt_Log_name,'file')),
    
    disp(['cant find file : ' opt_Log_name  ])
    error
else
    load  (opt_Log_name)
end


if notDefined('SunGrid');SunGrid=0;end
if notDefined('proclus');proclus=0;end
PDGridFit_done=false;

%%

fNum=ceil(length(opt.wh)/opt.jumpindex);
sgename=opt.SGE;
tic
while PDGridFit_done~=true
    % List all the files that have been created from the call to the
    % grid
    list=ls(opt.dirname);
    % Check if all the files have been made.  If they are, then collect
    % all the nodes and move on.
    if length(regexp(list, '.mat'))>=fNum,
        PDGridFit_done=true;
        % Once we have collected all the nodes we delete the sge outpust
        eval(['!rm -f ~/sgeoutput/*' sgename '*'])
    else
        % cheack if ther are jobs on the sun grid que list
        qStatCommand    = [' qstat | grep -i  job_' sgename(1:6)];
        [status result] = system(qStatCommand);
        tt=toc;
        if (isempty(result) && tt>60)
            % cheack if 1 min pass and there are no job waiting to be finish (and we don't have all the jobdone)
            %then we will need to re run it.
            
            RunSelectedJob=true;
            
            mrQ_fitM0boxesCall(opt_Log_name,SunGrid,proclus,RunSelectedJob)
           
            
        end
        
    end
end
