function mrQ_fitB1LR_Call(opt_logname,SunGrid,RunSelectedJob,GridOutputDir,clobber)
% function mrQ_fitB1LR_Call(opt_logname,SunGrid,RunSelectedJob,clobber)
%
%This function loads the opt structure that contains all the fit
% information and sends it to the computer grid if one is defined (SunGrid
% and/or proclus). If no grid is defined, the function will send the fit
% information to the local computer solver (SGE is faster, but it works
% fine without it as well). Depending on whether SunGrid is called, the
% function will either perform the voxel-wise fits using SunGrid or simply
% alone.
%
%  ~INPUTS~
%     opt_logname: Directory to the .mat file containing the appropriate
%                  parameters, as created in mrQ_PD_LRB1SPGR_GridParams.m
%         SunGrid: Whether to use the SunGrid; default is no.
%  RunSelectedJob: Default is false. If true, it will look for the missing
%                  jobs and then run them.
%         clobber: Overwrite existing data and reprocess. [Default = false]
%
%
% See also: mrQ_Gridcheck.m
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%

%% I. Load and definitions
load(opt_logname);
dirname=opt.dirname;
sgename=opt.SGE;
jumpindex=opt.jumpindex ;
fullID=sgename(isstrprop(sgename, 'digit'));
id=str2double(fullID(1:8));

if notDefined('clobber')
    clobber =false;
end

if clobber && (exist(dirname,'dir'))
    % In the case we start over and there are old fits,
    % so we will delete them
    eval(['! rm -r ' dirname]);
end

if (~exist(dirname,'dir')),
    mkdir(dirname);
end
if notDefined('GridOutputDir')
    GridOutputDir=pwd;
end
%% II.  Perform the gain fits
% II-a: With SunGrid

% Perform the fits for each box using the SunGrid Engine
if SunGrid==1; %if there is a grid
    
    if notDefined('RunSelectedJob')
        
        % Check to see if there is an existing SGE job that can be
        % restarted . if not, run all of them.
        
        % should we control the sgeoutput clear it before we start?
        %    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        %    sgerun2('mrQ_B1_LRFit(opt,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
        
        for jobindex=1:ceil(opt.N_Vox2Fit/jumpindex)
            jobname=1000*str2double(fullID(1:3))+jobindex;
            command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_B1_LRFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
            [stat,res]= system(command);
            
            if ~mod(jobindex,50)
                fprintf('%g jobs out of %g have been submitted   \n',jobindex,ceil(opt.N_Vox2Fit/jumpindex))
            end
        end
        
        %         sgerun('mrQ_B1_LRFit(opt_logname,jumpindex,jobindex);',sgename,1,1:ceil(opt.N_Vox2Fit/jumpindex),[],[],5000);
    else
        %         there is an existing job that can be restarted, ask user if to
        %         finish it ('yes') or restart the jobs from the beginning ('no').
        an1 = input( 'Unfinished SGE run found: Would you like to try and finish the existing SGE run? Press 1 if yes. To start over, press 0 ');
        if an1==1
            % Run only the selected jobs
            MissingFileNumber=mrQ_multiFit_WhoIsMissing( dirname,opt.N_Vox2Fit,jumpindex); % the job to run
            for kk=1:length(MissingFileNumber)
                jobindex=MissingFileNumber(kk);
                jobname=1000*str2double(fullID(1:3))+jobindex;
                command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_B1_LRFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                [stat,res]= system(command);
                
                if ~mod(kk,50)
                    fprintf('%g jobs out of %g have been submitted    \n',kk,length(MissingFileNumber));
                end
            end
        elseif an1==0
            for jobindex=1:ceil(opt.N_Vox2Fit/jumpindex)
                jobname=1000*str2double(fullID(1:3))+jobindex;
                command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_B1_LRFit(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                [stat,res]= system(command);
                
                if ~mod(jobindex,50)
                    fprintf('%g jobs out of %g have been submitted   \n',jobindex,ceil(opt.N_Vox2Fit/jumpindex))
                end
            end
        end
        %            sgerun('mrQ_B1_LRFit(opt_logname,jumpindex,jobindex);',sgename,1,MissingFileNumber,[],[],5000);
        
        % (We can also add an interactive call.
        % The code for that is commented below.)
        
    end
    
    
    % II-b: Without SunGrid
    
else %if there is no grid
    
    % Without a grid, this will take very long.
    disp('No parallel computation grid is used to fit PD. Using the local machine instead.');
    disp('      This may take a long time!!!');
    
    % In this case, the jumpindex is the number of voxels to fit (no jumps, only one job)
    jumpindex=   opt.N_Vox2Fit;
    opt.jumpindex=jumpindex;
    save(opt.logname,'opt');
    %     mrQ_B1_LRFit(opt_logname,jumpindex,1);
    mrQ_B1_LRFit(id,jumpindex,1);
    
    
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


