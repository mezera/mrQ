function mrQ_fitB1LR_Call(opt_logname,SunGrid,RunSelectedJob)
% this function load the opt straction that have all the fit information
% and send it in to the computer grid if one is defined (SunGrid and/or
% proclass).  If not it will send it to local compoter solver ( sge is
% faster but it is working fine with out it as well
%




load (opt_logname);
dirname=opt.dirname;
sgename=opt.SGE;
jumpindex=opt.jumpindex ;



%%   Perform the gain fits
% Perform the fits for each box using the Sun Grid Engine
if SunGrid==1;
    
     if notDefined('RunSelectedJob')
    
    % Check to see if there is an existing SGE job that can be
    % restarted. If not start the job, if yes prompt the user.
    if (~exist(dirname,'dir')),
        mkdir(dirname);
    end
    
    % should we control the sgeoutput clear it before we start?
    %    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        %    sgerun2('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
        
            sgerun('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
     else
         % run only the selected jobs this 
                      MissingFileNumber=mrQ_multiFit_HowIsMissing( dirname,N_Vox2Fit,jumpindex); % the job to  run
                   sgerun('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,MissingFileNumber,[],[],5000);
                   % we can also had an interactive call. the code for that
                   % is commented below.
     end
          
    
else
    % with out grid call that will take very long
    disp(  'No parallre computation grid is used to fit PD. Using the local machin instaed , this may take  long time !!!');
    % in this case the jumpindex is the number of voxel to fit (no jumps , only one job)
    jumpindex=   opt.N_Vox2Fit;
    opt.jumpindex=jumpindex;
    
    mrQ_B1_LRFit(opt,jumpindex,1);
        save(opt.logname,'opt');

    
    
end






%%

        
%         if notDefined('RunSelectedJob')
%             
%             % Prompt the user
%             inputstr = sprintf('An existing SGE run was found. \n Would you like to try and finish the exsist SGE run?');
%             RunSelectedJob = questdlg( inputstr,'mrQ_B1_LRFit','Yes','No','Yes' );
%             if strcmpi(RunSelectedJob,'yes'), RunSelectedJob = true; end
%             if strcmpi(RunSelectedJob,'no'),  RunSelectedJob = false; end
%         end
%         
%         % User opted to try to finish the started SGE run
%         if RunSelectedJob==true
%             MissingFileNumber=mrQ_multiFitControl( dirname,N_Vox2Fit,jumpindex)
%             
%             %% should we clean the SGE output Dir ?
%             if length(find(reval)) > 0
%                 eval(['!rm -f ~/sgeoutput/*' sgename '*'])
%             end
%         end
%          
%         
% %                     for kk=1:length(reval)
% %                         sgerun2('mrQ_B1_LRFit(opt,jumpindex,jobindex);',[sgename num2str(kk)],1,MissingFileNumber(kk),[],[],5000);
% %                     end
%                     sgerun('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,MissingFileNumber,[],[],5000);
%                 end
%             end
%             
%             % User opted to restart the existing SGE run
%         elseif RunSelectedJob==false,
%             t = pwd;
%             cd (opt.outDir)
%             eval(['!rm -rf ' dirname]);
%             cd (t);
%             eval(['!rm -f ~/sgeoutput/*' sgename '*'])
%             mkdir(dirname);
%             if proclass==1
%                 sgerun2('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
%             else
%                 sgerun('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
%                 
%             end
%         else
%             error('User cancelled');
%         end
%         
%         end
    
        
