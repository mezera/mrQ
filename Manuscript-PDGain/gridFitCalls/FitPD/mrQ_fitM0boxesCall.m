function mrQ_fitM0boxesCall(opt_logname,SunGrid,proclass,RunSelectedJob)
% this function load the opt straction that have all the fit information
% and send it in to the computer grid if one is defined (SunGrid and/or
% proclass).  If not it will send it to local compoter solver
%




load (opt_logname);
dirname=opt.dirname;
sgename=opt.SGE;
jumpindex=opt.jumpindex ;



%%   Perform the gain fits
% Perform the fits for each box using the Sun Grid Engine
if SunGrid==1;
    
    % Check to see if there is an existing SGE job that can be
    % restarted. If not start the job, if yes prompt the user.
    if (~exist(dirname,'dir')),
        mkdir(dirname);
        eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        if proclass==1
            sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
        else
            sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
            
        end
    else
        
        if notDefined('RunSelectedJob')
            
            % Prompt the user
            inputstr = sprintf('An existing SGE run was found. \n Would you like to try and finish the exsist SGE run?');
            RunSelectedJob = questdlg( inputstr,'mrQ_fitM0boxesCall','Yes','No','Yes' );
            if strcmpi(RunSelectedJob,'yes'), RunSelectedJob = true; end
            if strcmpi(RunSelectedJob,'no'),  RunSelectedJob = false; end
        end
        
        % User opted to try to finish the started SGE run
        if RunSelectedJob==true
            reval = [];
            list  = ls(dirname);
            ch    = 1:jumpindex:length(opt.wh);
            k     = 0;
            
            for ii=1:length(ch),
                ex=['_' num2str(ch(ii)) '_'];
                if length(regexp(list, ex))==0,
                    k=k+1;
                    reval(k)=(ii);
                end
            end
            
            if length(find(reval)) > 0
                eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                if proclass==1
                    for kk=1:length(reval)
                        sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',[sgename num2str(kk)],1,reval(kk),[],[],5000);
                    end
                else
                    sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,reval,[],[],5000);
                end
            end
            
            % User opted to restart the existing SGE run
        elseif RunSelectedJob==false,
            t = pwd;
            cd (opt.outDir)
            eval(['!rm -rf ' dirname]);
            cd (t);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            mkdir(dirname);
            if proclass==1
                sgerun2('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
            else
                sgerun('mrQ_CoilPD_gridFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex),[],[],5000);
                
            end
        else
            error('User cancelled');
        end
        
    end
    
else
    % with out grid call that will take very long
    disp(  'No parallre computation grid is used to fit PD. Using the local machin instaed , this may take very long time !!!');
    jumpindex=   length(opt.wh);
    opt.jumpindex=jumpindex;
    
    mrQ_CoilPD_gridFit(opt,jumpindex,1);
    save(opt.logname,'opt');
end

end