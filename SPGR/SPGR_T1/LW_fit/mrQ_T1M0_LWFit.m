function [T1_WL, T1_L,PD_WL, PD_L] = mrQ_T1M0_LWFit(s,brainMask,tr,flipAngles,Gain,B1,outDir,xform,SGE,savenow,sb,proclass)
%
% [T1 PD resNorm] = mrQ_T1M0_LWFit(s,brainMask,tr,flipAngles,M0,t1,Gain,B1,outDir,xform,SGE,savenow,sb,proclass)
%
% Perform weighted least square fitting of T1 and PD according to:
% Linear least-squares method for unbiased estimation of T1 from SPGR signals. Chang LC, Koay CG, Basser PJ, Pierpaoli C. Magn Reson Med. 2008 Aug;60(2):496-501.
%
%
% INPUTS:
%       s           - contains aligned data
%       brainMask   - Tissue mask delineating the brain region
%       tr          - TR taken from the S2 structure of aligned data
%       flipAngles  - Array of flipAngles for each scan.
%       M0          - MAP
%       t1          - [t1 data]
%       Gain        -
%       B1          -
%       outDir      - Ouput directory where the resulting nifti files will
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
%
%
% See Also:
%       mrQfit_T1M0_ver2.m
%
% (c) Stanford University, VISTA Lab
%
% Author: Aviv Mezer and Nikola Stikov date 24.02.2014 
% Copyright, Stanford University 2014
% 

%% Check inputs

if (~exist('sb','var')|| isempty(sb)),
    sb='UN';
end

sgename=[sb '_3dT1PD'];

if (~exist('SGE','var')|| isempty(SGE)),
    SGE=0;
end
if (~exist('proclass','var')|| isempty(proclass))
    proclass=0;
end

if (~exist('savenow','var')|| isempty(savenow)),
    savenow=0;
end

%% Set options for optimization procedure
a=version('-date');
if str2num(a(end-3:end))>=2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'off','Tolx',1e-12);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'off','Tolx',1e-12);%'TolF',1e-12
    
end% we put all the relevant data in a structure call op.t thiss will make it  easyer to send it between the computer in the grid
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

%% Perform the optimization (optionally using the SGE)

% USE THE SGE

clear brainMask tmp Res M0 options
if SGE==1;
    
    if (~exist([outDir '/tmpSGWL'],'dir')), mkdir([outDir '/tmpSGWL']);
        % the result form the grid will be saved in a tmporery directory
        if proclass==1
            sgerun2('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        else
            sgerun('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        end
    else
        an1 = input( 'Unfinished SGE run found: Would you like to try and finish the existing sge run? Press 1 if yes. To start over press 0 ');
        
        % Continue existing SGE run from where we left it last time
        % we find the fit that are missing
        if an1==1
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
                if proclass==1
                    a=num2str(ceil(rand(1)*10));
                    %sgerun2('mrQ_fitT1PD_SGE(opt,2000,jobindex);',[sgename a],1,reval)
                    for kk=1:length(reval)
                        sgerun2('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',[sgename num2str(kk)],1,reval(kk)); % we run the missing oupput again
                    end
                    
                else
                    sgerun('mrQ_fitT1PD_SGE(opt,8000,jobindex);',sgename,1,reval);
                end
            end
            list=ls(opt.outDir);
            
            % Restart the SGE processing from the beginning
        elseif an1==0
            t=pwd;
            cd (outDir)
            !rm -rf tmpSGWL
            cd (t);
            
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            mkdir([outDir '/tmpSGWL']);
            if proclass==1
                sgerun2('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
            else
                sgerun('mrQ_fitT1PDLW_SGE(opt,8000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
            end
        else
            error;
        end
    end
    
    %% build the data that was fit by the SGE to a T1 nd M0 maps
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
            
            
            qStatCommand    = [' qstat | grep -i  job_' sgename(1:6)];
            [status result] = system(qStatCommand);
            tt=toc;
            if (isempty(result) && tt>60)
                % then the are no jobs running we will need to re run it.
                
                %we will rerun only the one we need
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
                    if proclass==1
                        
                        % sgerun2('mrQ_fitT1PD_SGE(opt,2000,jobindex);',[sgename 'redo'],1,reval); % we run the missing oupput again
                        
                        for kk=1:length(reval)
                            sgerun2('mrQ_fitT1PDLW_SGE(opt,2000,jobindex);',[sgename num2str(kk)],1,reval(kk)); % we run the missing oupput again
                        end
                        
                    else
                        sgerun('mrQ_fitT1PDLW_SGE(opt,2000,jobindex);',sgename,1,reval);
                    end
                end
                
            else
                %  keep waiting
            end
            
            
            
            
            
        end
        
       
        
        
    end
    
    % NO SGE
    %using the local computer to fit T1 and the sunGrid
else
    
    fprintf('\n fit the T1 map localy, may be slow. SunGrid use can be much faster             \n');
    
    if (~exist([outDir '/tmpSGWL'],'dir')),
        mkdir([outDir '/tmpSGWL']);
        jobindex=1:ceil(length(opt.wh)/jumpindex);
    else
        jobindex=[];
        list=ls(opt.outDir);
        ch= 1:jumpindex:length(opt.wh) ;
        k=0;
        for ii=1:length(ch),
            
            ex=['_' num2str(ch(ii)) '_'];
            if length(regexp(list, ex))==0,
                k=k+1;
                jobindex(k)=(ii);
            end
        end
    end
    
    
    if ~isempty(jobindex)
        for i=jobindex
            mrQ_fitT1PDLW_SGE(opt,jumpindex,i);
        end
    end
    
    %Build the  T1 and M0 maps
    fNum=ceil(length(opt.wh)/jumpindex);
    % List all the files that have been created from the call to the
    % grid
    
    list=ls(opt.outDir);
    % Check if all the files have been made.  If they are, then collect
    % all the nodes and move on.
    
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


%% Save out results

if savenow==1
    dtiWriteNiftiWrapper(single(T1_WL), xform, fullfile(outDir,['T1_WL_last.nii.gz']));
    dtiWriteNiftiWrapper(single(PD_WL), xform, fullfile(outDir,['PD_WL_last.nii.gz']));
    dtiWriteNiftiWrapper(single(T1_L), xform, fullfile(outDir,['T1_L_last.nii.gz']));
    dtiWriteNiftiWrapper(single(PD_L), xform, fullfile(outDir,['PD_L_last.nii.gz']));
end