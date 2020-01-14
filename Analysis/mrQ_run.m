function mrQ_run(dir,outDir,inputData_spgr,inputData_seir,B1file, varArgIn)
% mrQ_run(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
%
%    example= mrQ_run_Ver2(dir,outdir,[],[],[],{'lsq',1})
%
%
%    INPUT:
%
%                dir:   Directory where the nifti from NIMS are located.
%             outDir:   Directory to which the output will be saved.
%                           (default: pwd/mrQ)
%     inputData_spgr:   The SPGR data
%     inputData_seir:   The SEIR data
%             B1file:   If empty (default), the function will calculate a
%                          B1 file from the data SEIR and SPGR data.
%                          Alternatively, the B1 inhomogeneity map can be
%                          provided as a NIfTI. The file has to be
%                          registered to the SPGR reference file. Note that
%                          if no reference file is provided, AC-PC
%                          alignment will be defined below inside the
%                          function mrQ_initSPGR_ver2.m.
%           varargin:   every parameter that you would like to be
%                          different from default can be given here. This
%                          input will be expected to come as a cell array
%                          and in pairs, according to the options given in
%                          the mrQ_set function. For example, if a
%                          reference file is given, the input should be
%                          {'ref', refFile} where refFile is the path to a
%                          reference image (nii.gz). Another example is the
%                          choice of using sungrid, in which case the input
%                          should be : {'sungrid',useSunGrid} where
%                          useSunGrid is 1/0 depending on whether or the
%                          user would like to use sungris (1) or not (0).
%                          the default is 0. One could also ask for both
%                          option: {'ref', refFile,'sungrid',useSunGrid}
%   example= mrQ_run_Ver2(dir,outdir,[],[],[],{'lsq',1})
%   OUTPUT:
%
%       This function creates and saves the mrQ strucure to the subject's
%       directory. New directories will be created, including directories
%       for data and quantitative fits. Images will be registered to each
%       other.
%            *  SEIR-EPI T1 will be computed (low resolution).
%            *  SPGR T1, M0, B1 maps, and a synthetic T1-weighted image
%                   will be computed.
%            *  T1-weighted and quantitative T1 images will be combined to
%                   segment the brain tissue.
%`           *  PD and coil gain will be fitted from the M0 image.
%            *  Biophysical models will be applied to calculate VIP and
%                   SIR maps.
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016
%
%

%% I. Create the initial structure

if notDefined('outDir')
    outDir = fullfile(dir,'mrQ');
end

% Creates the name of the output directory
if ~exist(outDir,'dir'); mkdir(outDir); end

% Creates the mrQ structure
mrQ = mrQ_Create(dir,[],outDir);

%     Create a file containing mrQpath, named after its 'ID' (taken from
%     its tempname). This allows for an easy use of SunGrid.
mrQ_createIDfile(mrQ);

% Set other parameters, such as SUNGRID and fieldstrength

if ~notDefined('varArgIn')
    if ~isempty(varArgIn)
        for ii = 1:2:numel(varArgIn)-1
            % Check to make sure that the argument is formatted properly
            mrQ = mrQ_Set(mrQ, varArgIn{ii}, varArgIn{ii+1});
        end
    end
end

%% II. Arrange the SPGR
% A specific arrange function for nimsfs, nifti, or using input for user

if ~isfield(mrQ,'ArrangeSPGR_Date');
    
    if ~notDefined('inputData_spgr')
        mrQ = mrQ_arrangeSPGR_nimsfs(mrQ,inputData_spgr);
    else
        mrQ = mrQ_arrangeSPGR_nimsfs(mrQ);
        
    end
else
    fprintf('Data was already arranged on %s \n',mrQ.ArrangeSPGR_Date)
end


%% IV. Initiate and align SPGR
%  parameters for aligning SPGR

%load(name);

if isfield(mrQ,'SPGR_init_done');
else
    mrQ.SPGR_init_done=0;
end

if     mrQ.SPGR_init_done==0
    
    % Keeps track of the variables we use.
    % For details, look inside the function.
    [mrQ.InitSPGR]=mrQ_initSPGR(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
    mrQ.SPGR_init_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n  init SPGR - done!           \n');
else
    fprintf(' \n Loading init SPGR data            \n');
    
end

%%  V. Fit SPGR PD

if ~isfield(mrQ,'SPGR_LinearT1fit_done');
    
    mrQ.SPGR_LinearT1fit_done=0;
end

% clobber is implemented inside (we can add this to the inputs)
if (mrQ.SPGR_LinearT1fit_done==0);
    
    [mrQ.LinFit]=mrQfit_T1M0_Lin(mrQ,mrQ.InitSPGR.spgr_initDir);
    
    mrQ.SPGR_LinearT1fit_done=1;
    
    save(mrQ.name,'mrQ');
    
    fprintf('\n Fit linear T1 SPGR  - done!              \n');
else
    fprintf('\n Loading linearly fitted SPGR T1                \n');
    
end

%% III. Perform SEIR fit
if ~notDefined('B1file')
    if exist(B1file,'file')
        mrQ.B1FileName=B1file;
        fprintf('Using the B1 map:  %s \n',B1file);
    else
        error('Can not find the B1 map:  %s \n',B1file);
    end
end

if  ~isfield(mrQ,'B1FileName')
    % Checks if B1 was defined by the user.
    % If not, we will use the SEIR data to map it.
    
    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    %%
    %
    %%
    if (mrQ.SEIR_done==0);
        
        % ARRANGE SEIR data
        if ~isfield(mrQ,'ArrangeSEIR_Date')
            
            if ~notDefined('inputData_seir')
                mrQ = mrQ_arrangeSEIR_nimsfs(mrQ,inputData_seir);
            else
                mrQ = mrQ_arrangeSEIR_nimsfs(mrQ);
            end
        else
            fprintf('Data was already arranged on %s \n',mrQ.ArrangeSEIR_Date)
        end
        
        %% fitiing SEIR T1 map and register to SPGE
        
        
        
        mrQ=mrQ_SEIR(mrQ);
        
    else
        fprintf('\n Loading previously fitted SEIR data \n');
    end
    
    
    %% VII. Fit & Build B1
    
    % Checks if B1 was defined by the user.
    
    if ~isfield(mrQ,'B1_done');
        mrQ.B1_done=0;
    end
    
    if ( mrQ.B1_done==0)
        
        
        mrQ=mrQ_B1_LR(mrQ);
        
        mrQ.B1_done=1;
        save(mrQ.name,'mrQ');
        fprintf('\n Building B1 - done!       \n');
        
    else
        fprintf(['Using existing B1 map file '   mrQ.B1FileName        '  \n']);
        
    end
end

%% VIII. T1M0 fit with B1

if ~isfield(mrQ,'SPGR_T1fit_done');
    mrQ.SPGR_T1fit_done=0;
end

if ( mrQ.SPGR_T1fit_done==0)
    
    mrQ=mrQ_Path_checks(mrQ);
    
    mrQ=mrQ_T1M0_Fit(mrQ);
    mrQ.SPGR_T1fit_done=true;
    save(mrQ.name,'mrQ');
    
    fprintf('\n fit T1 SPGR  - done!              \n');
    
else
    
    fprintf('\n Using previously fitted SPGR T1                \n');
    
end

%%  Segmentation needed for PD fit
% Prefer to PD fit
% 1. Get a segmentation (need freesurfer output)
% 2. Get CSF
% 3. Make a M0 file for the coils

%% IX. Create the synthetic T1-weighted images and save them to disk

if ~isfield(mrQ,'synthesis')
    mrQ.synthesis=0;
end
if mrQ.synthesis==0
    
    [mrQ.SegInfo.T1wSynthesis_MOT1,mrQ.SegInfo.T1wSynthesis_T1] =mrQ_T1wSynthesis1(mrQ,[],[],mrQ.LinFit.HeadMask);
    
    mrQ.synthesis=1;
    save(mrQ.name, 'mrQ')
    
    fprintf('\n Synthesis of T1  - done!              \n');
else
    fprintf('\n Using previously synthesized T1              \n');
end


%% X. Segmentation and CSF

if ~isfield(mrQ,'segmentation')
    mrQ.segmentation=0;
end

if mrQ.segmentation==0
    
    %     default- FSL segmentation
    if (mrQ.runfreesurfer==0 && ~isfield(mrQ,'freesurfer'))
        % Segment the T1w by FSL (step 1) and get the tissue mask (CSF-WM-GM) (step 2)
        %         mrQ=mrQ_segmentT1w2tissue(mrQ);
        mrQ=mrQ_Seg_kmeans(mrQ,mrQ.BrainMask);
        mrQ.segmentation=1;
        
        %      run FreeSurfer: it is slow and needs extra definitions.
    elseif (mrQ.runfreesurfer==1)
        mrQ=mrQ_Complitfreesurfer(mrQ);
        mrQ.segmentation=1;
        
        %      use an uploaded freesurfer nii.gz
    elseif   isfield(mrQ,'freesurfer');
        [mrQ.SegInfo]=mrQ_CSF(mrQ.spgr_initDir,mrQ.freesurfer,[],mrQ.AnalysisInfo);
        mrQ.segmentation=1;
    end
    
    save(mrQ.name,'mrQ');
    fprintf('\n Segmentation and CSF  - done!              \n');
else
    fprintf('\n Using previously segmented data              \n');
    
    
end

%% XI. Fitting PD from M0

if ~isfield(mrQ,'PDdone')
    mrQ.PDdone=0;
end
if mrQ.PDdone==0
    
    mrQ=mrQ_M0_ToPD(mrQ);
    mrQ.PDdone=1;
    save(mrQ.name, 'mrQ')
    fprintf('\n Calculation of PD from M0  - done!              \n');
else
    fprintf('\n Using previously calculated PD              \n');
end


%% XII. Calculate VIP, TV,  SIR and synthetic T1w

if ~isfield(mrQ,'VIP_WF_done')
    mrQ.VIP_WF_done=0;
end

if (mrQ.VIP_WF_done==0)
    fprintf('\n Calculate VIP, TV and SIR form T1 and WF maps               \n');
    
    [mrQ] = mrQ_WF(mrQ);
    
    
    [mrQ] = mrQ_VIP(mrQ);
    
    mrQ.VIP_WF_done=1;
    save(mrQ.name,'mrQ');
    
    fprintf('\n Calculation of VIP, MTV and SIR  - done!              \n');
    %
    % XIII. Create a series of synthetic T1w images
    
    [mrQ.T1w_file,mrQ.T1w_file1,mrQ.T1w_file2,mrQ.T1w_file3] =mrQ_T1wSynthesis2(mrQ);
    save(mrQ.name,'mrQ');
    fprintf('\n Calculation synthetic T1w- done!              \n');
    
else
    fprintf('\n Using previously calculated VIP, MTV, SIR and synthetic T1w             \n');
end




%% XIV. Organize the OutPut directory
if ~isfield(mrQ,'AnalysisDone')
    mrQ.AnalysisDone=0;
end

if (mrQ.AnalysisDone==0)
    mrQ=mrQ_arrangeOutPutDir(mrQ);
    mrQ.AnalysisDone=1;
    mrQ.AnalysisDoneDate=date;
end

mrQ_plotSummary(mrQ)
mrQ_deleteIDfile(mrQ);% delete the temporary ID file stored in mrQ/sge_subjects

%save
save(mrQ.name,'mrQ');
fprintf('\n done !!! \n')
