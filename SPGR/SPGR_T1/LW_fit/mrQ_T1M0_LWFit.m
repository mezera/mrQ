function [T1_WL, T1_L,PD_WL, PD_L] = mrQ_T1M0_LWFit(s,brainMask,tr, ...
    flipAngles,Gain,B1,outDir,xform,mrQ,GridOutputDir,savenow)
%
%[T1_WL, T1_L,PD_WL, PD_L] = mrQ_T1M0_LWFit(s,brainMask,tr, ...
%                            flipAngles,Gain,B1,outDir,xform,mrQ,savenow)
%
%Performs weighted least square fitting of T1 and PD according to: Linear
% least-squares method for unbiased estimation of T1 from SPGR signals.
% Chang LC, Koay CG, Basser PJ, Pierpaoli C. Magn Reson Med. 2008
% Aug;60(2):496-501.
%
%
% INPUTS:
%       s           - contains aligned data
%       brainMask   - Tissue mask delineating the brain region
%       tr          - TR taken from the S2 structure of aligned data
%       flipAngles  - Array of flipAngles for each scan.
%       M0          - M0 map
%       t1          - T1 data
%       Gain        -
%       B1          -
%       outDir      - Ouput directory where the resulting NIfTI files will
%                     be saved.
%       xform       - Transform
%       SGE         - Option to run using SGE [default = 0]
%       savenow     - Saves the outputs to disk [default = 0]
%       sub         - Subject name for SGE call
%
%
% OUTPUTS:
%       B1
%       resNorm
%       PD
%
%
% WEB RESOURCES
%
% See also: mrQfit_T1M0_ver2.m
%
% (c) Stanford University, VISTA Lab
%
% Authors: Aviv Mezer and Nikola Stikov, 24.02.2014
% Copyright, Stanford University, 2014
%

%% I. Check inputs

% if (~exist('sb','var')|| isempty(sb)),
%     sb='UN';
% end
%
% sgename=[sb '_3dT1PD'];

% if (~exist('SGE','var')|| isempty(SGE)),
%     SGE=0;
% end
% if (~exist('proclass','var')|| isempty(proclass))
%     proclass=0;
% end
if (~exist('savenow','var')|| isempty(savenow)),
    savenow=0;
end
if isfield(mrQ,'SunGrid')
    SGE=mrQ.SunGrid;
else
    SGE=0;
end

if isfield(mrQ,'sub')
    sb=mrQ.sub;
else
    sb='UN';
end
sgename=[sb '_3dT1PD'];

fullID=sb(isstrprop(sb, 'digit'));
id=str2double(fullID(1:8));

if notDefined('GridOutputDir')
    GridOutputDir=pwd;
end
%% II. Set options for optimization procedure
% Adjust for different versions of matlab
a=version('-date');
if str2num(a(end-3:end))>=2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'off','Tolx',1e-12);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'off','Tolx',1e-12);%'TolF',1e-12
end
% We put all the relevant data in a structure called "opt".
% This will make it easier to send it between the computers in the grid

sz=size(brainMask);

for i=1:length(s)
    
    tmp=s(i).imData(brainMask);
    opt.s(:,i)=double(tmp);
end

opt.flipAngles = double(flipAngles);
opt.tr = double(tr);

opt.wh   = find(brainMask);
opt.Gain = double(Gain(brainMask));
opt.B1   = double(B1(brainMask));

opt.outDir = [outDir '/tmpSGWL'];
opt.name   = '/T1PDlsqVx';
jumpindex=8000;


% Save a logfile with all the options used during processing:
logname = [outDir '/fitT1LW.mat'];
opt.logname=logname;

% Save an information file we can load afterwards, if needed.
save(opt.logname,'opt');
%added this save of mrQ
mrQ.LWoptname=logname;
save(mrQ.name,'mrQ');

%% III. Perform the optimization (optional: use the SunGrid)
%% III-a. USE THE SGE

clear brainMask tmp Res M0 options
if SGE==1;
    
    if (~exist([outDir '/tmpSGWL'],'dir')), mkdir([outDir '/tmpSGWL']);
        % the result form the grid will be saved in a temporary directory
        for jobindex=1:ceil(length(opt.wh)/jumpindex)
            jobname=1000*str2double(fullID(1:3))+jobindex;
            command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_fitT1PDLW_SGE(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
            [stat,res]=system(command);
            if ~mod(jobindex,100)
                fprintf('%g jobs out of %g have been submitted         \n',jobindex,ceil(length(opt.wh)/jumpindex))
            end
            
        end
        %         if proclass==1
        %             sgerun2('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        %         else
        %             sgerun('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        %         end
        
    else   % if there's an existing SGE job, ask user if to restart this job or start from scratch
        an1 = input( 'Unfinished SGE run found: Would you like to try and finish the existing SGE run? Press 1 if yes. To start over, press 0 ');
        
        % Continue existing SGE run from where we left it last time
        % we find the fits that are missing
        if an1==1
            MissingFileNumber=mrQ_multiFit_WhoIsMissing(opt.outDir,length(opt.wh),jumpindex);
            if ~isempty(MissingFileNumber)
                
                for kk=1:length(MissingFileNumber)
                    jobindex=MissingFileNumber(kk);
                    jobname=num2str(100*str2double(fullID(1:3)));
                    command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_fitT1PDLW_SGE(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                    [stat,res]= system(command);
                    if ~mod(kk,100)
                        fprintf('%g jobs out of %g have been submitted       \n',kk,length(MissingFileNumber));
                    end
                end
            end
            
            
        elseif an1==0   % Restart the SGE processing from the beginning
            t=pwd;
            cd (outDir)
            !rm -rf tmpSGWL
            cd (t);
            
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            mkdir([outDir '/tmpSGWL']);
            for jobindex=1:ceil(length(opt.wh)/jumpindex)
                jobname=1000*str2double(fullID(1:3))+jobindex;
                command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_fitT1PDLW_SGE(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                [stat,res]= system(command);
                if ~mod(jobindex,100)
                    fprintf('%g jobs out of %g have been submitted           \n',jobindex,ceil(length(opt.wh)/jumpindex));
                end
                
            end
            %             if proclass==1
            %                 sgerun2('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
            %             else
            %                 sgerun('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
            %             end
        else
            error('user cancelled');
        end
    end
    
    % Build the data that was fit by the SGE to a T1 and M0 maps
    
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
            
            % Loop over the nodes and collect the output
            for i=1:fNum
                st=1 +(i-1)*jumpindex;
                ed=st+jumpindex-1;
                
                if ed>length(opt.wh), ed=length(opt.wh);end
                
                name=[opt.outDir '/' opt.name '_' num2str(st) '_' num2str(ed) '.mat'];
                load (name);
                t11(st:ed)=res(1,:);
                pd1(st:ed)=res(2,:);
                t12(st:ed)=res(3,:);
                pd2(st:ed)=res(4,:);
                
            end
            % Once we have collected all the nodes we delete the temporary
            t=pwd;
            cd (outDir)
            !rm -r tmpSGWL
            cd (t);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            
        else
            qStatCommand    = [' qstat | grep -i  job_' fullID(1:3)];
            [status result] = system(qStatCommand);
            tt=toc;
            if (isempty(result) && tt>60)
                % then there are no jobs running we will need to re run it.
                MissingFileNumber=mrQ_multiFit_WhoIsMissing(opt.outDir,length(opt.wh),jumpindex);
                if length(find(MissingFileNumber))>0
                    % clean the sge output dir and run the missing fit
                    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                    for kk=1:length(MissingFileNumber)
                        jobindex=MissingFileNumber(kk);
                        jobname=num2str(100*str2double(fullID(1:3)));
                        command=sprintf('qsub -cwd -j y -b y -N job_%g "matlab -nodisplay -r ''mrQ_fitT1PDLW_SGE(%f,%g,%g); exit'' >log"', jobname, id,jumpindex,jobindex);
                        [stat,res]= system(command);
                        if ~mod(kk,100)
                            fprintf('%g jobs out of %g have been submitted       \n',kk,length(MissingFileNumber));
                        end
                    end
                end
                
            else
                %  keep waiting
            end
            
        end
        
    end
    
    %% III-b: NO SunGrid
    % using the local computer to fit T1 and PD
else
    
    fprintf('\n Fitting the T1 map locally, may be slow. SunGrid use can be much faster             \n');
    
    % Check for matlab version and for parallel computing toolbox
    MyVer = ver; % check matlab version   
    has_PCTbox = any(strcmp(cellstr(char(MyVer.Name)), 'Parallel Computing Toolbox')); % check for PCTbox
    MyVer_ed=MyVer.Release; % identify release version
    MyVer_year= sscanf(MyVer_ed,'%*[^0123456789]%d'); % identify release year
    MyVer_AorB= sscanf(MyVer_ed,'%*[^ab]%c'); % identify version a or b

    
    if (~exist([outDir '/tmpSGWL'],'dir')),
        mkdir([outDir '/tmpSGWL']);
        MissingFileNumber=1:ceil(length(opt.wh)/jumpindex); % all of the files are 'missing'
    else
        MissingFileNumber= mrQ_multiFit_WhoIsMissing( opt.outDir,length(opt.wh),jumpindex);
        %         jobindex=[];
        %         list=ls(opt.outDir);
        %         ch= 1:jumpindex:length(opt.wh) ;
        %         k=0;
        %         for ii=1:length(ch),
        %
        %             ex=['_' num2str(ch(ii)) '_'];
        %             if length(regexp(list, ex))==0,
        %                 k=k+1;
        %                 jobindex(k)=(ii);
        %             end
        %         end
    end
 
  % Parallel Processing
  % Using parallel processing can reduce the runtime of the LW fit,
  % from ~60 minutes to ~20 minutes
  
  if has_PCTbox == 0 %no PCTbox, and thus no parfor
      if ~isempty(MissingFileNumber)
          for kk=1:length(MissingFileNumber) %regular for-loop
              jobindex=MissingFileNumber(kk);
              mrQ_fitT1PDLW_SGE(id,jumpindex,jobindex);
          end
      end
  else %PCTbox exists, and so does parfor
      if MyVer_year<2013 || MyVer_year==2013 && MyVer_AorB=='a' % check if matlab is 2013a or earlier
          if ~isempty(MissingFileNumber)
              % Find number of available workers
              myworkers=findResource;
              myworkers=myworkers.ClusterSize;
              
              matlabpool('open', myworkers); %opens available pools (matlab 2013a and earlier)
              parfor kk=1:length(MissingFileNumber) %parallel for-loop
                  
                  jobindex=MissingFileNumber(kk);
                  mrQ_fitT1PDLW_SGE(id,jumpindex,jobindex);
              end
              matlabpool close; %close pools
          end
      elseif MyVer_year>2013 || MyVer_year==2013 && MyVer_AorB=='b'% check if matlab is 2013b or later
          if ~isempty(MissingFileNumber)
              tic
              gcp(); %opens available pools (matlab 2013b and later)
              parfor kk=1:length(MissingFileNumber) %parallel for-loop
                  jobindex=MissingFileNumber(kk);
                  mrQ_fitT1PDLW_SGE(id,jumpindex,jobindex);
              end
              delete(gcp); %close pools
              toc
          end
      end
  end
    
    
    
%     if ~isempty(MissingFileNumber)
%         for kk=1:length(MissingFileNumber)
%             jobindex=MissingFileNumber(kk);
%             mrQ_fitT1PDLW_SGE(id,jumpindex,jobindex);
%         end
%     end
    
    %Build the  T1 and M0 maps
    fNum=ceil(length(opt.wh)/jumpindex);
    
    % Loop over the nodes and collect the output
    for i=1:fNum
        st=1 +(i-1)*jumpindex;
        ed=st+jumpindex-1;
        
        if ed>length(opt.wh), ed=length(opt.wh);end
        
        name=[opt.outDir '/' opt.name '_' num2str(st) '_' num2str(ed) '.mat'];
        load (name);
        t11(st:ed)=res(1,:);
        pd1(st:ed)=res(2,:);
        t12(st:ed)=res(3,:);
        pd2(st:ed)=res(4,:);
    end
    % Once we have collected all the nodes we delete the temporary
    t=pwd;
    cd (outDir)
    !rm -r tmpSGWL
    cd (t);
    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
    
end

T1_WL = zeros(sz);
PD_WL = T1_WL; T1_L=T1_WL;PD_L = T1_WL;
T1_WL(opt.wh) = t11(:)./1000;
PD_WL(opt.wh) = pd1(:);
T1_L(opt.wh) = t12(:)./1000;
PD_L(opt.wh) = pd2(:);

%% IV. Save out results

if savenow==1
    dtiWriteNiftiWrapper(single(T1_WL), xform, fullfile(outDir,['T1_WL_last.nii.gz']));
    dtiWriteNiftiWrapper(single(PD_WL), xform, fullfile(outDir,['PD_WL_last.nii.gz']));
    dtiWriteNiftiWrapper(single(T1_L), xform, fullfile(outDir,['T1_L_last.nii.gz']));
    dtiWriteNiftiWrapper(single(PD_L), xform, fullfile(outDir,['PD_L_last.nii.gz']));
end

if SGE
        jobname=fullID(1:3);
    filesPath=[GridOutputDir,'/job_',num2str(jobname),'*'];
    delCommand=sprintf('rm %s', filesPath);
    [status, result]=system(delCommand);
end