function mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
%  mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
%
% This is an improved version of: mrQ_runNIMS(dir,Callproclus,refFile,outDir)
%
%    INPUT:
%
%                dir:   Directory where the nifti from NIMS are located.
%             outDir:   Directory to which the output will be saved. 
%                           (default: pwd/mrQ)
%         useSUNGRID:   Whether to use SUNGRID cluster computing
%            refFile:   The path to a reference image (nii.gz)
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
%
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
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%  2015
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

% Set other parameters, such as SUNGRID and fieldstrength

%            mrQ = mrQ_Set(mrQ,'sub',num2str(ii));

if notDefined('useSUNGRID')
    mrQ = mrQ_Set(mrQ,'sungrid',false);
else
    mrQ = mrQ_Set(mrQ,'sungrid',useSUNGRID);
end

%     Create a file containing mrQpath, named after its 'ID' (taken from
%     its tempname). This allows for an easy use of SunGrid.
mrQ_createIDfile(mrQ);

    
%             mrQ = mrQ_Set(mrQ,'sungrid',1);
mrQ = mrQ_Set(mrQ,'fieldstrength',3);

if ~notDefined('refFile')
    mrQ = mrQ_Set(mrQ,'ref',refFile);
else
    % New input to automatically perform ac-pc alignment
    mrQ = mrQ_Set(mrQ,'autoacpc',1);
end

%% II. Arrange the data
% A specific arrange function for nimsfs, nifti, or using input for user

if ~isfield(mrQ,'Arrange_Date');
    
    if (~notDefined('inputData_spgr') &&  ~notDefined('inputData_seir'))
        mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir);
    else
        mrQ = mrQ_arrangeData_nimsfs(mrQ);
        
    end
else
    fprintf('Data was already arranged on %s \n',mrQ.Arrange_Date)
end
%% III. Perform SEIR fit

if notDefined('B1file')
    % Checks if B1 was defined by the user.
    % If not, we will use the SEIR data to map it.
    
    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    
    if (mrQ.SEIR_done==0);
        
        % Keeps track of the variables we use.  
        % For details, see inside the function.
        [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR_ver2(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag);
        
        [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,[],0,mrQ);
        
        mrQ.SEIR_done=1;
        save(mrQ.name,'mrQ');
        fprintf('Fit SEIR  - done! \n');
        
    else
        fprintf('\n Loading previously fitted SEIR data \n');
        
    end
    
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
    [~, ~, ~,~,~, mrQ]=mrQ_initSPGR_ver2(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
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
    
    [mrQ]=mrQfit_T1M0_Lin(mrQ);
    
    mrQ.SPGR_LinearT1fit_done=1;
    
    save(mrQ.name,'mrQ');
    
    fprintf('\n Fit linear T1 SPGR  - done!              \n');
else
    fprintf('\n Loading linearly fitted SPGR T1                \n');
    
end

%%  VI. Register high-resolution EPI image to low-resolution aligned T1 image

%mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)

if ~isfield(mrQ,'SPGR_EPI_align_done');
    mrQ.SPGR_EPI_align_done=0;
end

if ( mrQ.SPGR_EPI_align_done==0)
    
    mrQ.spgr2epiAlignFile=fullfile(mrQ.spgr_initDir,'SEIRepiSPGRAlign_best_RB.mat');
    [mrQ.Ants_Info]=mrQ_NLANTS_warp_SPGR2EPI_RB(mrQ.SEIR_epi_T1file, mrQ.T1_LFit_HM, mrQ.SPGR_niiFile_FA, mrQ.spgr_initDir, mrQ.spgr2epiAlignFile, mrQ.AligndSPGR);
    
    mrQ.SPGR_EPI_align_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n Alignment of EPI to T1  - done!              \n');
else
    fprintf(['\n Using alignment of EPI to T1, calculated on '    mrQ.Ants_Info. spgr2epi_Align_date           '\n']);
    
end

%% VII. Build B1

if notDefined('B1file')
    % Checks if B1 was defined by the user.
    
    if ~isfield(mrQ,'B1Build_done');
        mrQ.B1Build_done=0;
    end
    
    if ( mrQ.B1Build_done==0)
        
        
        mrQ=mrQ_B1_LR(mrQ);
        
        mrQ.B1Build_done=1;
        save(mrQ.name,'mrQ');
        fprintf('\n Building B1 - done!       \n');
        
    else
        fprintf(['Using existing B1 map file '   mrQ.B1FileName        '  \n']);
        
    end
else
    mrQ.B1FileName=B1file;
end

%% VIII. T1M0 fit with B1

if ~isfield(mrQ,'SPGR_T1fit_done');
    mrQ.SPGR_T1fit_done=0;
end

if ( mrQ.SPGR_T1fit_done==0)
    
    
    mrQ=mrQ_T1M0_Fit(mrQ);
    mrQ.SPGR_T1fit_done=true;
    save(mrQ.name,'mrQ');
    
    fprintf('\n fit T1 SPGR  - done!              \n');
    
else
    
    fprintf('\n using previously fitted SPGR T1                \n');
    
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
    
    [mrQ.SegInfo.T1wSynthesis_MOT1,mrQ.SegInfo.T1wSynthesis_T1] =mrQ_T1wSynthesis1(mrQ,[],[],mrQ.HeadMask);
    
    mrQ.synthesis=1;
    save(mrQ.name, 'mrQ')
    
    fprintf('\n Synthesis of T1  - done!              \n');
else
    fprintf('\n using previously synthesized T1              \n');
end


%% X. Segmentation and CSF

if ~isfield(mrQ,'segmentation');
    mrQ.segmentation=0;
end

if mrQ.segmentation==0;
    
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
    fprintf('\n segmentation and CSF  - done!              \n');
else
    fprintf('\n using previously segmented data              \n');
    
    
end

%% XI. Fitting PD from M0

if ~isfield(mrQ,'PDdone')
    mrQ.PDdone=0;
end
if mrQ.PDdone==0
    
   mrQ=mrQ_M0_ToPD(mrQ);
    
    save(mrQ.name, 'mrQ')
    mrQ.PDdone=1;
    fprintf('\n Calculation of PD from M0  - done!              \n');
else
    fprintf('\n Using previously calculated PD              \n');
end

%% XII. Calculate VIP, TV, and SIR

if ~isfield(mrQ,'SPGR_PDBuild_done')
    mrQ.SPGR_PDBuild_done=0;
end

if (mrQ.SPGR_PDBuild_done==0)
    fprintf('\n Calculate VIP, TV and SIR form T1 and WF maps               \n');
    
    [mrQ.AnalysisInfo] = mrQ_VIP(mrQ);
    save(mrQ.name,'mrQ');
    mrQ.SPGR_PDBuild_done=1;
    fprintf('\n Calculation of VIP, TV and SIR  - done!              \n');
    
else 
     fprintf('\n Using previously calculated VIP, TV and SIR              \n');
end


%%  XIII. Create a series of synthetic T1w images

[mrQ.T1w_file,mrQ.T1w_file1] =mrQ_T1wSynthesis1(mrQ);

%% XIV. Organize the OutPut directory
mrQ=mrQ_arrangeOutPutDir(mrQ);
mrQ_deleteIDfile(mrQ);% delete the temporary ID file stored in mrQ/sge_subjects

%done
mrQ.AnalysisDone=1;
mrQ.AnalysisDoneDate=date;
%save
save(mrQ.name,'mrQ');
fprintf('\n done !!! \n')
