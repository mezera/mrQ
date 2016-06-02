function mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ)
% mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ)
%
% This function will try to register the SPGR T1 to a fitted T1 SEIR using ANTS.
% We will make sure that the registration is good by the gradient between
% the two images.
% in cases the registration is off up to a threshold (mrQ.QuantAntsThresh).
% Note that this threshold might be different between scanner and dipend on image %quality and epi artifacts.  WE find 0.65 to be good on  GE Discovery MR750 and 0.8 for Siemens Skyra
% we will try different starting point.
% To do so we will register the SEIr image each time to a different
% inversion time image.
% Additionally, we will try to register the SPGR to an ACPC space or to native space.
% Last we will crop the SPGR image and remove dark area that can effect the
% registration.
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016

if isfield(mrQ,'QuantAntsThresh');
else
    mrQ.QuantAntsThresh=0.65;
end

Ants_Info.QuantAnts2High(1:length(mrQ.inputdata_seir.IT),1:3)=inf;

for jj=1:3
    if jj==1
        SPGR_T1_file=mrQ.LinFit.T1_LFit_HM;
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
    elseif jj==2
        
        if isnan(mrQ.refIm) % if so we are allready in native space. no need
        else
            [mrQ.InitSPGR_NtSp]=mrQ_initSPGR(mrQ.SPGR,nan,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
            spgr_initDir=mrQ.InitSPGR_NtSp.spgr_initDir;
            
            [mrQ.LinFit_NtSp]=mrQfit_T1M0_Lin(mrQ,spgr_initDir);
            
            SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit_HM;
            SPGRFieldName='InitSPGR_NtSp';
            
            mrQ.B1Raw2ACPC=true;
            save(mrQ.name,'mrQ');
            
        end
    elseif jj==3
        if isfield(mrQ,'NBF_Tresh');
        else
            mrQ.NBF_Tresh=0.75;
        end
        
        SPGR_T1_file=mrQ.LinFit.T1_LFit_HM;
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        
        t1=readFileNifti(SPGR_T1_file);
        t1.data(:,:,find(mrQ.LinFit.FracNonBrainT1<mrQ.NBF_Tresh))=0;
        t1fileHMCliped = fullfile(spgr_initDir,['T1_LFit_HMCliped.nii.gz']);
        dtiWriteNiftiWrapper(single(t1.data), t1.qto_xyz, t1fileHMCliped);
        mrQ.LinFit.t1fileHMCliped=t1fileHMCliped;
        SPGR_T1_file=t1fileHMCliped;
        mrQ.B1Raw2ACPC=false;
        mrQ.B1ClipT1=true;
        
        mrQ.Ants_Info=Ants_Info;
        save(mrQ.name,'mrQ');
        
        
    end
    
    if ~isfield (mrQ,'SEIRfitDone'); mrQ.SEIRfitDone(1:length(mrQ.inputdata_seir.IT))=0;end
    
    
    
    for ii=1: length(mrQ.inputdata_seir.IT) %% in each iteration, align the SEIR images with a different image as the target image
        % Keeps track of the variables we use.
        % For details, see inside the function.
        
        Dir_SEIR_Align2Im = fullfile(mrQ.outDir, ['SEIR_Align2Im_' num2str(ii)]);
        
        SEIR_fit_dir = fullfile(Dir_SEIR_Align2Im, 'fitT1_GS');
        
        
        if mrQ.SEIRfitDone(ii)==0;
            
            mkdir(SEIR_fit_dir);
            
            [~, ~, ~, SEIRsaveData]=mrQ_initSEIR(mrQ,Dir_SEIR_Align2Im,mrQ.alignFlag,ii);
            
            %SEIR_epi_AlignIm=mrQ.inputdata_seir.name{ii};  %keep the image that is the target image in the current iteration.
            
            [mrQ.SEIRfits{ii}]=mrQ_fitSEIR_T1(Dir_SEIR_Align2Im,[],0);
            mrQ.SEIRfits{ii}.Dir_SEIR_Align2Im=Dir_SEIR_Align2Im;
            mrQ.SEIRfits{ii}.SEIR_fit_dir=SEIR_fit_dir;
            
            mrQ.SEIRfitDone(ii)=1;
                    mrQ.Ants_Info=Ants_Info;

            save(mrQ.name,'mrQ');
            
            
        end
        %%   Register high-resolution EPI image to low-resolution aligned T1 image
        
        %mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)
        %
        % if ~isfield(mrQ,'SPGR_EPI_align_done');
        %     mrQ.SPGR_EPI_align_done=0;
        % end
        %
        % if ( mrQ.SPGR_EPI_align_done==0)
        
        % create a temp folder for the ...
        
        AntsPath = fullfile(SEIR_fit_dir,['ANTS_' num2str(jj)]);
        if ~exist(AntsPath,'dir'), mkdir(AntsPath); end
        
        if isfield(mrQ,'antsmet')
            %             if mrQ.antsmet=='1'
            %                 % fitting based on the brain only- no head
            %                 [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR{ii}.SEIR_epi_T1file,mrQ.T1_LFit,mrQ.SEIR_epi_Maskfile,AntsPath,{mrQ.T1_LFit_HM});
            %             elseif mrQ.antsmet=='2'
            %                 % fitting without a mask on the epi
            %                 [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR{ii}.SEIR_epi_T1file,mrQ.T1_LFit_HM,[],AntsPath,{mrQ.T1_LFit_HM});
            %             end
        else
            
            [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIRfits{ii}.SEIR_epi_T1file,SPGR_T1_file,mrQ.SEIRfits{ii}.SEIR_epi_Maskfile,AntsPath,{SPGR_T1_file});
        end
        
        
        MovingScaleConstat=1000;% to compare SPGR T1 in sec to SEIR T1 in ms.
        % check the registration:
        T1_spgr_epi = T1_spgr_epi{1};
        [Ants_Info.QuantAntsScore]=mrQ_QuantAnts(mrQ.SEIRfits{ii}.SEIR_epi_T1file,T1_spgr_epi,MovingScaleConstat);
        % We will keep the files only if the are the current best Ants registration
        if min(Ants_Info.QuantAnts2High(:)) > Ants_Info.QuantAntsScore
            
            
            
            Ants_Info.WARP_SPGR_EPI = WARP_SPGR_EPI;
            
            Ants_Info.T1_spgr_epi=T1_spgr_epi;
            
            % Saving the epi alignment parameters
            
            
            
            Ants_Info.SEIR_SPGR_Curent_AlignNums=[ii,jj];
            Ants_Info.SEIR_SPGR_Curent_AlignDirs={SEIR_fit_dir, spgr_initDir};
            Ants_Info.SPGRFieldName=SPGRFieldName;
        end
        
        
        if Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh;
                        Ants_Info.QuantAnts2High(ii,jj)=0;

            %saving the correct SEIR files
            
            mrQ.SEIR_epi_T1file = mrQ.SEIRfits{ii}.SEIR_epi_T1file;
            
            mrQ.SEIR_epi_resnormfile = mrQ.SEIRfits{ii}.SEIR_epi_resnormfile;
            
            mrQ.SEIR_epi_fitFile = mrQ.SEIRfits{ii}.SEIR_epi_fitFile;
            
            mrQ.SEIR_epi_M0file = mrQ.SEIRfits{ii}.SEIR_epi_M0file;
            
            mrQ.SEIR_epi_Maskfile =mrQ.SEIRfits{ii}.SEIR_epi_Maskfile;
            
            % deleting temp folders
            %         cmd = ['rm -Rf ', AntsTmpPath];
            %         system(cmd);
            %        we can cosider removing all the extra epi fits;
            %         cmd = ['rm -Rf ', SEIRFitTmpPath];
            %         system(cmd);
            %
                    mrQ.Ants_Info=Ants_Info;

                        save(mrQ.name,'mrQ');

            break
        else
            % if the registration quality did not pass the threshold, we will
            % keep the values of the registration quality and move on to the
            % next iteration.
            
            %         Ants_Info.WARP_SPGR_EPI = best_WARP_SPGR_EPI;
            %         Ants_Info.T1_spgr_epi = best_T1_spgr_epi;
            Ants_Info.QuantAnts2High(ii,jj)=Ants_Info.QuantAntsScore;
            %  mrQ.SEIR_epi_AlignImToHigh{ii}=SEIR_epi_AlignIm;
                    mrQ.Ants_Info=Ants_Info;

            save(mrQ.name,'mrQ');
        end
    end
    
    if Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh;        break; end
    
    
    
end