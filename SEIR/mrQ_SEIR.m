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
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016


mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ);
% % % % % if ANTS was unsuccessfull in registring SPGR and EPI, TRY fsl
% % % % if mrQ.Ants_Info.QuantAntsScore > mrQ.QuantAntsThresh;   
% % % %     mrQ=mrQ_Call_FSLAlign_forSEIR_SPGR(mrQ);
% % % % end
% IF FSL was also unsuccessfull, break and look at the data
if mrQ.Ants_Info.QuantAntsScore > mrQ.QuantAntsThresh  && mrQ.testR1_BM==1  
    error('we can not trust the EPI-SPGR registration \nPlease manually check the registration between \n %s and \n %s \n If it is ok, manually change mrQ.SEIR_done to be =1, save mrQ, and rerun.', mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi)
elseif mrQ.Ants_Info.QuantAntsScore > mrQ.QuantAntsThresh  && mrQ.testR1_BM==0
    error('we can not trust the EPI-SPGR registration \nPlease manually check the registration between \n %s and \n %s \n If it is ok, manually change mrQ.SEIR_done to be =1, save mrQ, and rerun. \n \n You can also test if the skull strip generated from R1: \n %s \n is better compred to the skull strip used in this run: \n %s \n If it is better, change mrQ.testR1_BM to be =1, save mrQ, and rerun', mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi,mrQ.SEIR_epi_T1bmfile_fromR1,mrQ.SEIR_epi_T1bmfile)
end

%% 
mrQ.SEIR_done=1;

% mrQ.SPGR_EPI_align_done=1;

save(mrQ.name,'mrQ');
fprintf('\n Alignment of EPI to T1  - done!              \n');




