function mrQ=mrQ_Complitfreesurfer(mrQ);
% comleeting the tfreesurfer segmentation and use it to make tissue masks
if isfield(mrQ,'freesurfer')
    % we got freesurfer 
    fprintf('\n load freesurfer seqmwntation file                \n ');
    
    
else %we need to get the freesufer segmentation here or to run it now
    c=clock;
    fprintf(['\n' date  ' ' num2str(c(4)) ':'  num2str(c(5))  ' working to get  the freesurfer segmentation; this can take time               \n'  ])
    subjID=['freesufer_' mrQ.AnalysisInfo.sub];
    
    if isfield(mrQ.AnalysisInfo,'freeSufer_subdir')  %freesurfer of the syntetic T1w image (was was start  inside mrQfit_T1M0_ver2.m)
        
        location=mrQ.AnalysisInfo.freeSufer_subdir; %freesufer write to here
        lstfilemripath=[location '/mri/wmparc.mgz']; %the last file that need to be made byfreeSufer
    
    elseif   isfield(mrQ,'freeSuferRefIm_subdir') % in the case freesurfer was done on the refIm that was provided
        
        location=mrQ.freeSuferRefIm_subdir;         %freesufer write to here
        lstfilemripath=[location '/mri/wmparc.mgz'];  %the last file that need to be made byfreeSufer
    
        if exist(lstfilemripath,'file') % got the last file so we  can move on with freesufer
            
            [mrQ.freesurfer mrQ.freesurfer1  ]= mrQfinalized_autosegment(subjID,mrQ.AnalysisInfo.T1wSynthesis,1);
    
        else
            % we don't have the file let check if it's still working or the
            % freesurfer process just fail along te way.
            
            wait_on=1;
             fprintf('\n freeSufer files are missing wating try to wait for it to finish (maybe long)              \n')
            
            while wait_on
                %cheak if anything was wrrithing in the last 5 hours or we will
                %run it again
            
                [status, result] =system([' find ' location ' -type f -mmin -300']);
                
                if      exist(lstfilemripath,'file') % got the last file so we  can move on
                
                    [mrQ.freesurfer mrQ.freesurfer1  ]=mrQfinalized_autosegment(subjID,mrQ.AnalysisInfo.T1wSynthesis,1);
                    wait_on=0;
                
                elseif  notDefined('result')
                    %we wait no file was writen and the last needed file is not there so somthing is wrong let run freesufer again
                    wait_on=0;
                    cmd=(['! rm -r ' location]); %clear the directory
                    eval(cmd);
                    [mrQ.freesurfer mrQ.freesurfer1 ] =  mrQfinalized_autosegment(subjID,mrQ.AnalysisInfo.T1wSynthesis);
                    
                end
            end
        end
    else
        
        %it freesufer was never startedso it have to be down now - this
        %may take about 24 hours
       fprintf('\n  freeSufer files are missing. starting segmentation now this maybe long ~24h               \n')
      
        [mrQ.freesurfer mrQ.freesurfer1  ]   =   mrQfinalized_autosegment(subjID,mrQ.AnalysisInfo.T1wSynthesis);
    end
end
save(mrQ.name,'mrQ');

% 2. CSF
if isfield(mrQ,'FITCSF');
else
    mrQ.FITCSF=0;
end

if (mrQ.FITCSF==0);
    fprintf('finding the CSF               ');
    
    % use thefreesurfer  segmentation
    [mrQ.AnalysisInfo]=mrQ_CSF(mrQ.spgr_initDir,mrQ.freesurfer);
    mrQ.FITCSF=1;
    save(mrQ.name,'mrQ');
else
    fprintf('load the CSF               ');
end

