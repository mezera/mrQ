function mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ)
% mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ)
%
% This function will try to register the SPGR T1 to a fitted T1 SEIR using
% ANTS. We will make sure that the registration is good by the gradient
% between the two images. in cases the registration is off up to a
% threshold (mrQ.QuantAntsThresh). Note that this threshold might be
% different between scanners and depends on image quality and epi
% artifacts.  We find 0.65 to be good on GE Discovery MR750 and 0.8 for
% Siemens Skyra. We will try different starting point. To do so we will
% register the SEIR image each time to a different inversion time image.
% Additionally, we will try to register the SPGR to an ACPC space or to
% native space. Last we will crop the SPGR image and remove dark area that
% can effect the registration.
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016

if ~isfield(mrQ,'QuantAntsThresh')
    mrQ.QuantAntsThresh=0.65;
end
if ~isfield(mrQ,'ants_bm')
    mrQ.ants_bm=0;
end

if ~isfield(mrQ,'testR1_BM')
    mrQ.testR1_BM=0;
end

% In some cases we get better BrainMask for the SEIR map with R1.
%if this flag is on we will try registration with the BrainMask from R1.
testR1_BM=mrQ.testR1_BM;

tryR1mask=0; % for other iterations don't use BrainMask from R1.
use_TtSp=0;
if testR1_BM==1  % if testR1_BM is on there will be 3 trials
    iter_num=3;
else
    iter_num=2;  % normally there will be be 2 trials
end

iter_start =1;

if isfield(mrQ,'Ants_Info') % If you rerun with flag testR1_BM and the previous iterations allready exists
    if isfield(mrQ.Ants_Info,'QuantAnts2High') && testR1_BM==1
        if any(mrQ.Ants_Info.QuantAnts2High(:)>0 & mrQ.Ants_Info.QuantAnts2High(:)<inf)
            iter_start=3;
            Ants_Info.QuantAnts2High=mrQ.Ants_Info.QuantAnts2High;
        end
    end
end
Ants_Info.QuantAnts2High(1:length(mrQ.inputdata_seir.IT),iter_start:iter_num)=inf;

for jj=iter_start:iter_num
    if jj==1
        %         in the first round of trials we are trying to realign the SEIR
        %         images in different orders and refit T1.
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        useBM=mrQ.ants_bm;
        msg='refit SEIR-T1 map ';
        
    elseif jj==2
        %          for the second round of trials, we intitally tried to clip the T1
        %         map in the SPGR, but this didn't seem to work. So now we try to
        %         register with/without a bm - the opposite of what was done until
        %         now.
        
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        useBM=~(mrQ.ants_bm);
        msg='change the BM input ';
        
%     elseif jj==3
%         %         In the third round of trials, in case we initially used a
%         %         refernce image for the PGR or aligned them to acpc space, we now
%         %         keep them in native space which is generally closer to the
%         %         SEIR-EPI data
%         if isempty(mrQ.refIm) || ~isnan(mrQ.refIm) % if so we are already in native space. no need
%             
%             use_TtSp=1;
%             useBM=mrQ.ants_bm;
%             mrQ.tmpRef=mrQ.refIm;
%             mrQ.refIm=nan;
%             [mrQ.InitSPGR_NtSp]=mrQ_initSPGR(mrQ.SPGR,[],mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
%             mrQ.refIm= mrQ.tmpRef;
%             spgr_initDir=mrQ.InitSPGR_NtSp.spgr_initDir;
%             
%             [mrQ.LinFit_NtSp]=mrQfit_T1M0_Lin(mrQ,spgr_initDir);
%             
%             %  SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit_HM;
%             SPGRFieldName='InitSPGR_NtSp';
%             
%             mrQ.B1Raw2ACPC=true;
%             save(mrQ.name,'mrQ');
%             msg='realign SPGR in it native space and refit T1, ';
%         else
%             continue
%         end
    elseif jj==3
        %         In the forth round of trials, if the flag testR1_BM==1, we try to register with a bm
        %         from R1 and not M0 (somtimes this mask is better).
        spgr_initDir=mrQ.InitSPGR.spgr_initDir;
        SPGRFieldName='InitSPGR';
        tryR1mask=1;
        msg='refit SEIR-T1 map with brain mask from R1 ';
        useBM=1;
%     elseif jj==5
%         %         If the flag testR1_BM==1 and ants_bm==1, and we initially used a
%         %         refernce image for the PGR or aligned them to acpc space, we now
%         %         keep them in native space which is generally closer to the
%         %         SEIR-EPI data and use the Brain Mask from R1 for SEIR
%         if  isempty(mrQ.refIm) || ~isnan(mrQ.refIm) % if so we are already in native space. no need
%             useBM=1;
%             tryR1mask=1;
%             use_TtSp=1;
%             mrQ.tmpRef=mrQ.refIm;
%             mrQ.refIm=nan;
%             [mrQ.InitSPGR_NtSp]=mrQ_initSPGR(mrQ.SPGR,[],mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
%             mrQ.refIm= mrQ.tmpRef;
%             spgr_initDir=mrQ.InitSPGR_NtSp.spgr_initDir;
%             
%             [mrQ.LinFit_NtSp]=mrQfit_T1M0_Lin(mrQ,spgr_initDir);
%             
%             %         SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit_HM;
%             SPGRFieldName='InitSPGR_NtSp';
%             
%             mrQ.B1Raw2ACPC=true;
%             save(mrQ.name,'mrQ');
%             msg='realign SPGR in it native space and refit T1 with  brain mask from R1 for SEIR ';
%         else
%             continue
%         end
    end
    
    if ~isfield (mrQ,'SEIRfitDone'); mrQ.SEIRfitDone(1:length(mrQ.inputdata_seir.IT))=0;end
    
    if jj>1
        display(sprintf(['The SPGR - EPI registration is not good enough \nWe will ',msg,' and then try again.']));
    end
    
    for ii=1: length(mrQ.inputdata_seir.IT)
        % in each iteration, align the SEIR images with a different image as the target image
        % Keeps track of the variables we use.
        % For details, see inside the function.
        
        if jj==1 & ii>1
            display(sprintf(['The SPGR - EPI registration is not good enough \nWe will ',msg,' and then try again.']));
        elseif jj>1 & ii>1
            display(sprintf(['The SPGR - EPI registration is not good enough \nWe will try with another pair of maps.']));
        end
        
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
        
        % create a temp folder for the ...
        
        AntsPath = fullfile(SEIR_fit_dir,['ANTS_' num2str(jj)]);
        if ~exist(AntsPath,'dir'), mkdir(AntsPath); end
        
        
        if useBM==1    % Register based on brain only- skull stripped
            SPGR_T1_file=mrQ.LinFit.T1_LFit;
            SEIR_T1file=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile;
            
            %  [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIRfits{ii}.SEIR_epi_T1file,SPGR_T1_file,mrQ.SEIRfits{ii}.SEIR_epi_Maskfile,AntsPath,{SPGR_T1_file});
        elseif useBM==0   % Register without EPI mask - with skull
            SPGR_T1_file=mrQ.LinFit.T1_LFit_HM;
            SEIR_T1file = mrQ.SEIRfits{ii}.SEIR_epi_T1file;
        end    % [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIRfits{ii}.SEIR_epi_T1file,SPGR_T1_file,[],AntsPath,{SPGR_T1_file});
        
        if tryR1mask==1    % Register based on brain only- skull stripped using R1
            SEIR_T1file=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile_fromR1;
            SPGR_T1_file=mrQ.LinFit.T1_LFit;
        end
        
        if use_TtSp==1 && useBM==0
            SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit_HM;
        elseif use_TtSp==1 && useBM==1
            SPGR_T1_file=mrQ.LinFit_NtSp.T1_LFit;
        end
        
        [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(SEIR_T1file,SPGR_T1_file,[],AntsPath,{SPGR_T1_file});
        
        MovingScaleConstat=1000;% to compare SPGR T1 in sec to SEIR T1 in ms.
        % check the registration:
        T1_spgr_epi = T1_spgr_epi{1};
        [Ants_Info.QuantAntsScore]=mrQ_QuantAnts(SEIR_T1file,T1_spgr_epi,MovingScaleConstat);
        %         [Ants_Info.QuantAntsScore]=mrQ_QuantAnts(mrQ.SEIRfits{ii}.SEIR_epi_T1file,T1_spgr_epi,MovingScaleConstat);
        
        % We will keep the files only if the are the current best Ants registration
        if min(Ants_Info.QuantAnts2High(:)) > Ants_Info.QuantAntsScore
            
            
            Ants_Info.WARP_SPGR_EPI = WARP_SPGR_EPI;
            Ants_Info.T1_spgr_epi=T1_spgr_epi;
            
            % Saving the epi alignment parameters
            Ants_Info.SEIR_SPGR_Curent_AlignNums=[ii,jj];
            Ants_Info.SEIR_SPGR_Curent_AlignDirs={SEIR_fit_dir, spgr_initDir};
            Ants_Info.SPGRFieldName=SPGRFieldName;
            
            %saving the SEIR files
            mrQ.SEIR_epi_T1file = mrQ.SEIRfits{ii}.SEIR_epi_T1file;
            mrQ.SEIR_epi_resnormfile = mrQ.SEIRfits{ii}.SEIR_epi_resnormfile;
            mrQ.SEIR_epi_fitFile = mrQ.SEIRfits{ii}.SEIR_epi_fitFile;
            mrQ.SEIR_epi_M0file = mrQ.SEIRfits{ii}.SEIR_epi_M0file;
            mrQ.SEIR_epi_Maskfile =mrQ.SEIRfits{ii}.SEIR_epi_Maskfile;
            if tryR1mask==1
                mrQ.SEIR_epi_Maskfile=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile_fromR1;
            end

            mrQ.SEIR_epi_T1bmfile_fromR1=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile_fromR1;
            mrQ.SEIR_epi_T1bmfile=mrQ.SEIRfits{ii}.SEIR_epi_T1bmfile;
            Ants_Info.Use_BrainMask = useBM;
        end
        
        %         check if the current (best) registration is good enough
        if Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh
            Ants_Info.QuantAnts2High(ii,jj)=0;
            
            % % deleting temp folders
            %         cmd = ['rm -Rf ', AntsTmpPath];
            %         system(cmd);
            % %       we can cosider removing all the extra epi fits;
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
    
    if Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh
        break;
    end
    
end