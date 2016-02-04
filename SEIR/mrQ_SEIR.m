function mrQ=mrQ_SEIR(mrQ)

if isfield(mrQ,'AntSQantTresh');  %DOC%!!!!! add to mrQ.init and mrQ.set
else
    mrQ.AntSQantTresh=0.65;
end


for ii=1: length(mrQ.inputdata_seir.IT) %%DOC%%
    % Keeps track of the variables we use.
    % For details, see inside the function.
    [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag,ii);
    
    mrQ.SEIR_epi_AlignIm=mrQ.inputdata_seir.name{ii};  %DOC%%
    
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
    
    [mrQ.AntSQant]=mrQ_QuantAnts(mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi,MovingScaleConstat);
    
    if mrQ.AntSQant< mrQ.AntSQantTresh;
        break
    else
        mrQ.AntSQantToHigh(ii)=mrQ.AntSQant;
        mrQ.SEIR_epi_AlignImToHigh{ii}=mrQ.SEIR_epi_AlignIm;
    end
end


if mrQ.AntSQant> mrQ.AntSQantTresh;
    error('we can not trust the EPI SPGR registration please check WE NEED TO FIX THIS NOTE...')
end


mrQ.SEIR_done=1;

% mrQ.SPGR_EPI_align_done=1;

save(mrQ.name,'mrQ');
fprintf('\n Alignment of EPI to T1  - done!              \n');




