function mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
% mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir)
%mrQ_runNIMS(dir,Callproclus,refFile,outDir)
%dir - where the nifti from NIMS are.
% % % % %Callproclus use 1 when using proclus (stanfrod computing cluster)
%refFile a path to a reference image (nii.gz)
%outDir  where the output will be saved (defult is: pwd/mrQ)

%


%%
% Create the initial structure
if notDefined('outDir')
    outDir = fullfile(dir,'mrQ');
end %creates the name of the output directory

if ~exist(outDir,'dir'); mkdir(outDir); end

mrQ = mrQ_Create(dir,[],outDir); %creates the mrQ structure

% Set other parameters
%            mrQ = mrQ_Set(mrQ,'sub',num2str(ii));

if notDefined('useSUNGRID')
    mrQ = mrQ_Set(mrQ,'sungrid',false);
else
    mrQ = mrQ_Set(mrQ,'sungrid',useSUNGRID);
end

%             mrQ = mrQ_Set(mrQ,'sungrid',1);
mrQ = mrQ_Set(mrQ,'fieldstrength',3);

if ~notDefined('refFile')
    mrQ = mrQ_Set(mrQ,'ref',refFile);
    
else
    % New input to automatically acpc align
    mrQ = mrQ_Set(mrQ,'autoacpc',1);
end

%% arrange data
% Specific arrange function for nimsfs nifti or using input for user
if ~isfield(mrQ,'Arange_Date');
    
    if (~notDefined('inputData_spgr') &&  ~notDefined('inputData_seir'))
        mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir);
    else
        mrQ = mrQ_arrangeData_nimsfs(mrQ);
        
    end
else
    fprintf('data was already arranged at %s \n',mrQ.Arange_Date)
end
%% fit SEIR   

if notDefined('B1file') 
    % check if B1 was defined by the user. if not we will use the SEIR data
    % to map it.

    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    
    if (mrQ.SEIR_done==0);
        
        %keep track of the variable we use  for detail see inside the function
        [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR_ver2(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag);
        
        [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,[],0,mrQ);
        
        mrQ.SEIR_done=1;
        save(mrQ.name,'mrQ');
        fprintf('fit SEIR  - done!');
        
    else
        fprintf('\n  load fit SEIR data ! \n');
        
    end
    
end
    
    
%% intiate and  Align SPGR
%  param for Align SPGR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load(name);

if isfield(mrQ,'SPGR_init_done');
else
    mrQ.SPGR_init_done=0;
end

if     mrQ.SPGR_init_done==0
    
    %keeps track of the variables we use;  for details, look inside the function
    [~, ~, ~,~,~, mrQ]=mrQ_initSPGR_ver2(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
    mrQ.SPGR_init_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n  init SPGR - done!           \n');
else
    fprintf(' \n load init SPGR data            \n');
    
end


%%  Fit SPGR PD

if ~isfield(mrQ,'SPGR_LinearT1fit_done');

    mrQ.SPGR_LinearT1fit_done=0;
end

% clobber is implemented inside (we can add this to the inputs)
if (mrQ.SPGR_LinearT1fit_done==0);
   
      [mrQ]=mrQfit_T1M0_Lin(mrQ);
%     [mrQ.AnalysisInfo]=mrQfit_T1M0_ver2(mrQ);
    
    mrQ.SPGR_LinearT1fit_done=1;
    
    save(mrQ.name,'mrQ');
    
    fprintf('\n fit linear T1 SPGR  - done!              \n');
else
    fprintf('\n load linearly fitted SPGR T1                \n');
    
end

%%  register high res EPI image to low res aligned T1 image

%mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)

if ~isfield(mrQ,'SPGR_EPI_align_done');
    mrQ.SPGR_EPI_align_done=0;
end

if ( mrQ.SPGR_EPI_align_done==0)
    
    mrQ.spgr2epiAlignFile=fullfile(mrQ.spgr_initDir,'SEIRepiSPGRAlign_best_RB.mat');
    [mrQ.Ants_Info, Res]=mrQ_NLANTS_warp_SPGR2EPI_RB(mrQ.SEIR_epi_T1file, mrQ.T1_LFit_HM, mrQ.SPGR_niiFile_FA, mrQ.spgr_initDir, mrQ.spgr2epiAlignFile, mrQ.AligndSPGR);
    
    mrQ.SPGR_EPI_align_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n alignment of EPI to T1  - done!              \n');
else
    fprintf(['\n using alignment of EPI to T1, precalculated on '    mrQ.Ants_Info. spgr2epi_Align_date           '\n']);
    
end


%% B1
if ~isfield(mrQ,'B1Build_done');
    mrQ.B1Build_done=0;
end

if ( mrQ.B1Build_done==0)
    

 mrQ=mrQ_B1_LR(mrQ);

    else
        fprintf(['Using the  B1  map  file '   mrQ.B1FileName        '  \n');
        
    end
    

% 
%% T1M0 fit with B1
if ~isfield(mrQ,'SPGR_T1fit_done');
    mrQ.SPGR_T1fit_done=0;
end

if ( mrQ.SPGR_T1fit_done==0)
    
    
    mrQ=mrQ_T1M0_Fit(mrQ);
    mrQ.SPGR_T1fit_done=true;
    fprintf('\n fit  T1 SPGR  - done!              \n');
    
else
    
    fprintf('\n load  fitted SPGR T1                \n');
    
end

%%  segmentation needed for PD fit
%prefer to PD fit 1. get a segmentation (need freesurfer output) 2. get CSF; 3.make a M0 fies for the coils

%%  Create the synthetic T1 weighted images and save them to disk

% currently not working

if ~isfield(mrQ,'synthesis')
    mrQ.synthesis=0;
end
if mrQ.synthesis==0

    % [mrQ.SegInfo.T1wSynthesis,mrQ.SegInfo.T1wSynthesis1] =mrQ_T1wSynthesis1(mrQ,mrQ.WFfile,mrQ.T1file,mrQ.BrainMask);
[mrQ.SegInfo.T1wSynthesis_MOT1,mrQ.SegInfo.T1wSynthesis_T1] =mrQ_T1wSynthesis1(mrQ,mrQ.M0_LFit_HM,mrQ.T1_LFit_HM,mrQ.HeadMask,mrQ.spgr_initDir);

mrQ.synthesis=1;
save(mrQ.name, 'mrQ')

fprintf('\n Synthesis of T1  - done!              \n');
else 
    fprintf('\n using previously synthesized T1              \n');
end

% [mrQ.SegInfo.T1wSynthesis,mrQ.SegInfo.T1wSynthesis1, mrQ.SegInfo.maskSynthesis]=mrQ_T1wSynthesis(dataDir,B1file,outDir);
% save(infofile,'AnalysisInfo');

%%
%. Segmentaion and CSF
if ~isfield(mrQ,'segmentaion');
    mrQ.segmentaion=0;
end

if mrQ.segmentaion==0;
    
    %     default- fsl segmentation
    if (mrQ.runfreesurfer==0 && ~isfield(mrQ,'freesurfer'))
        % Segment the T1w by FSL (step 1) and get the tissue mask (CSF WM GM) (step 2)
%         mrQ=mrQ_segmentT1w2tissue(mrQ,BMfile,T1file,t1wfile,outDir,csffile,boxsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this doesn't work yet because the function is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calling non existent fields and files
        mrQ=mrQ_segmentT1w2tissue(mrQ,[],mrQ.SegInfo.T1wSynthesis_T1);
        mrQ.segmentaion=1;
        
        %      run FreeSurfer : it is slow and needs extra defintions.
    elseif (mrQ.runfreesurfer==1)
        mrQ=mrQ_Complitfreesurfer(mrQ);
        mrQ.segmentaion=1;
        
        %      use an uploaded freesurfer nii.zg
    elseif   isfield(mrQ,'freesurfer');
        [mrQ.SegInfo]=mrQ_CSF(mrQ.spgr_initDir,mrQ.freesurfer,[],mrQ.AnalysisInfo);
        mrQ.segmentaion=1;
        
    end
    save(mrQ.name,'mrQ');
    fprintf('\n segmentation and CSF  - done!              \n');
else
    fprintf('\n using previously segmented data              \n');
    
    
end

%% fitting PD from M0

M0file=WEneedAFunction;


% set parameters for the fit
[mrQ.opt_logname]=mrQ_PD_multicoil_RgXv_GridCall_Fittesting(mrQ.spgr_initDir,mrQ.sub,mrQ.PolyDeg,M0file,T1file,BMfile,[],[],[],Coilsinfo,Reg,Init,[],mrQ);


%%
%%
% calculate VIP TV and SIR

if (mrQ.SPGR_PDBuild_done==1)
    fprintf('\n calculate VIP TV SIR form T1 and WF maps               \n');
    
    [mrQ.AnalysisInfo] = mrQ_VIP(mrQ);
    save(mrQ.name,'mrQ');
end





%% Organize the OutPut  directory
mrQ=mrQ_arangeOutPutDir(mrQ);

%%  Create a series of synthetic T1w images

 [mrQ.T1w_file,mrQ.T1w_file1] =mrQ_T1wSynthesis1(mrQ);
%done
mrQ.AnalysisDone=1;
mrQ.AnalysisDoneDate=date;
%save
save(mrQ.name,'mrQ');
fprintf('\n done !!! \n')