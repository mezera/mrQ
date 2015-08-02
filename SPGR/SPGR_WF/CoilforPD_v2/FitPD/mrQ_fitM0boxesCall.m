function mrQ_fitM0boxesCall(opt_logname,SunGrid,RunSelectedJob)
% this function load the opt straction that have all the fit information
% and send it in to the computer grid if one is defined (SunGrid and/or
% proclass).  If not it will send it to local compoter solver
%


load (opt_logname);
dirname=opt.dirname;
sgename=opt.SGE;
jumpindex=opt.jumpindex ;

fullID=sgename(isstrprop(sgename, 'digit'));
id=str2double(fullID(1:8));

%%   Perform the gain fits
% Perform the fits for each box using the Sun Grid Engine
if SunGrid==1;
    jumpindex=5;
    
    % Check to see if there is an existing SGE job that can be
    % restarted. If not start the job, if yes prompt the user.
    if (~exist(dirname,'dir')),
        mkdir(dirname);
        eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        %         if proclass==1
        %             sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
        %         else
        %             sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
        %
        %         end
        for jobindex=1:ceil(length(opt.wh)/jumpindex)
            jobname=1000*str2double(fullID(1:3))+jobindex;
            command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_CoilPD_gridFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
            [stat,res]= system(command);
              
              if ~mod(jobindex,50)
                  fprintf('%g jobs out of %g have been submitted',jobindex,ceil(length(opt.wh)/jumpindex));
              end
        end
        
    else
        
        if notDefined('RunSelectedJob')
            
            % Prompt the user
            inputstr = sprintf('An existing SGE run was found. \n Would you like to try and finish the exist SGE run?');
            RunSelectedJob = questdlg( inputstr,'mrQ_fitM0boxesCall','Yes','No','Yes' );
            if strcmpi(RunSelectedJob,'yes'), RunSelectedJob = true; end
            if strcmpi(RunSelectedJob,'no'),  RunSelectedJob = false; end
        end
        
        % User opted to try to finish the started SGE run
        if RunSelectedJob==true
            reval = []; % will hold indices of chunks to re-evaluate
            list  = ls(dirname);
            ch    = 1:jumpindex:length(opt.wh); % Beginning indices of each chunk(?)
            fullID     = 0;
            
            for ii=1:length(ch), % loop over all chuncks, to see which needs re-evaluation
                ex=['_' num2str(ch(ii)) '_'];
                if length(regexp(list, ex))==0,
                    fullID=fullID+1;
                    reval(fullID)=(ii);
                end
            end
            
            if length(find(reval)) > 0
                eval(['!rm -f ~/sgeoutput/*' sgename '*'])
%                 if proclass==1
                    for kk=1:length(reval)
                        jobindex=reval(kk);
                        jobname=1000*str2double(fullID(1:3))+jobindex;
                        %   sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',[sgename num2str(kk)],1,reval(kk),                       [],[],5000);
                        command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_CoilPD_gridFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                        [stat,res]=system(command);
                        if ~mod(kk,50)
                            fprintf('%g jobs out of %g have been submitted',kk,length(reval));
                        end
                    end
%                 else
%                     sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,reval,[],[],5000);
%                 end
            end
            
            % User opted to restart the existing SGE run
        elseif RunSelectedJob==false,
            t = pwd;
            cd (opt.outDir)
            eval(['!rm -rf ' dirname]);
            cd (t);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            mkdir(dirname);
            for jobindex=1:ceil(length(opt.wh)/jumpindex)
                jobname=100*str2double(fullID(1:3))+jobindex;
                command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_CoilPD_gridFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                [stat,res]=system(command);
                if ~mod(jobindex,50)
                  fprintf('%g jobs out of %g have been submitted',jobindex,ceil(length(opt.wh)/jumpindex));
                end
            end
            %             if proclass==1
            %                 sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
            %             else
            %                 sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
            %
            %             end
        else
            error('User cancelled');
        end
        
    end
    
      % This loop checks if all the outputs have been saved and waits until
    % they are all done
    StopAndSave=0;
    fNum=ceil(length(opt.wh)/jumpindex);
    tic
    while StopAndSave==0
        % List all the files that have been created from the call to the
        % grid
        
        list=ls(opt.outDir);
        % Check if all the files have been made.  If they are, then collect
        % all the nodes and move on.
        if length(regexp(list, '.mat'))==fNum,
            StopAndSave=1;
        else
            jobname=fullID(1:3);
            qStatCommand    = [' qstat | grep -i  job_' jobname];
            [status result] = system(qStatCommand);
            tt=toc;
            if (isempty(result) && tt>60)
                % then the are no jobs running we will need to re run it.
                
                %we will rerun only the ones we need
                reval=[];
                list=ls(opt.outDir);
                ch=[1:jumpindex:length(opt.wh)];
                k=0;
                for ii=1:length(ch),
                    
                    ex=['_' num2str(ch(ii)) '_'];
                    if length(regexp(list, ex))==0,
                        k=k+1;
                        reval(k)=(ii);
                    end
                end
                
                if length(find(reval))>0
                    % clean the sge output dir and run the missing fit
                    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                    
                    for kk=1:length(reval)
                        jobindex=reval(kk);
                        jobname=1000*str2double(fullID(1:3))+jobindex;
                        %   sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',[sgename num2str(kk)],1,reval(kk),                       [],[],5000);
                        command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_CoilPD_gridFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                        [stat,res]=system(command);
                        if ~mod(kk,50)
                            fprintf('%g jobs out of %g have been submitted',kk,length(reval));
                        end
                   
                    end
                   
                end
                
            else
                %  keep waiting
            end
               
        end
    end
        
        
    
else
    % with out grid call that will take very long
    disp(  'No parallel computation grid is used to fit PD. Using the local machine instead , this may take a very long time !!!');
    if (~exist(dirname,'dir')),
        mkdir(dirname);
    end
    jumpindex=   length(opt.wh);
    opt.jumpindex=jumpindex;
    
    mrQ_CoilPD_gridFit(id,jumpindex,1);
    save(opt.logname,'opt');
end

end
