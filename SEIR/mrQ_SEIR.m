function mrQ=mrQ_SEIR(mrQ)
% we found the registration is often failing (causing a bias in the T1 fit
% in following steps). A possible solution is to re-align (and then re-fit)
% the SEIR, with a different IT image as the target of the registration. To
% monitor that, we register according to the first TI, and then check the
% registration using the funcction mrQ_QuantAnts. The check is not optimal
% and is using a somewhat arbitrary threshold. if the check doesn't pass
% the threshold - we reallign, fit and register spgr to the new fit - and
% then check again. if non of the trials pass our threshold, we ask the
% user to mannually check the registration, and if it is satisfying one can
% easilly bypass out checks using, for example, rerun mrQ with the
% additional input {'seir_done',1}: 
% mrQ_run(inputDir,outDir,[],[],[], {'seir_done',1}),

%%
if isfield(mrQ,'AntSQantTresh');  
else
    mrQ.AntSQantTresh=0.65;
end


for ii=1: length(mrQ.inputdata_seir.IT) %% in each iteration, align the SEIR images with a different image as the target image
    % Keeps track of the variables we use.
    % For details, see inside the function.
    [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag,ii);
    
    mrQ.SEIR_epi_AlignIm=mrQ.inputdata_seir.name{ii};  %keep the image that is the target image in the current iteration. 
    
    [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,[],0,mrQ);
    
    save(mrQ.name,'mrQ');
    fprintf('Fit SEIR  - done! \n');
    
    
    
    %%   Register high-resolution EPI image to low-resolution aligned T1 image
    
    %mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)
    %
    % if ~isfield(mrQ,'SPGR_EPI_align_done');
    %     mrQ.SPGR_EPI_align_done=0;
    % end
    %
    % if ( mrQ.SPGR_EPI_align_done==0)
    
    
    [mrQ.Ants_Info.WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR_epi_T1file,mrQ.T1_LFit_HM,mrQ.SEIR_epi_Maskfile,mrQ.spgr_initDir,{mrQ.T1_LFit_HM});
    mrQ.Ants_Info.T1_spgr_epi=T1_spgr_epi{1};
   
   MovingScaleConstat=1000;% to compare SPGR T1 is sec to SEIR T1 in ms.
    % check the registration:
    [mrQ.AntSQant]=mrQ_QuantAnts(mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi,MovingScaleConstat);
    
    if mrQ.AntSQant< mrQ.AntSQantTresh;
        break
    else
        % if the registration quality did not pass the threshold, we will
        % keep the values of the registration quality and move on to the
        % next iteration. 
        mrQ.AntSQantToHigh(ii)=mrQ.AntSQant;
        mrQ.SEIR_epi_AlignImToHigh{ii}=mrQ.SEIR_epi_AlignIm;
    end
end


if mrQ.AntSQant> mrQ.AntSQantTresh;
    error('we can not trust the EPI-SPGR registration \nPlease mannually check the registration between \n %s and \n %s \n If it is ok, mannually adjust the mrQ.mrQ.AntSQantTresh accordingly and rerun mrQ', mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi)
end


mrQ.SEIR_done=1;

% mrQ.SPGR_EPI_align_done=1;

save(mrQ.name,'mrQ');
fprintf('\n Alignment of EPI to T1  - done!              \n');




