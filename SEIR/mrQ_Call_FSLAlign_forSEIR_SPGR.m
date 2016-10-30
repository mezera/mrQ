
function  mrQ=mrQ_Call_FSLAlign_forSEIR_SPGR(mrQ)

% this function is used as a last resort for the registration of epi to
% spgr. It uses FSL using unix call, and iterates over the different T1
% fits to find the best registration. Her ewe first use linear
% registratuion to bring the maps close together, and then warp the map to
% get a better match.  This might increase the smoothing in the fit and
% should be recosidered (the 1-step warping didn't work for me)

%% more options for FSL registration are in the function:
% /home/shai.berman/Documents/Code/testing_pipelines/registration/fslFnirt/testFnirt.m
%% 
display(sprintf(['Ants was unable to register SPGR - EPI maps.\nWe will use FSL, and try the different pairs of maps.']));


SPGR_T1_file=mrQ.LinFit.T1_LFit;% mrQ.LinFit.T1_LFit_HM; maybe sing the skull stripped map is better
spgr_initDir=mrQ.InitSPGR.spgr_initDir;
SPGRFieldName='InitSPGR';
jj=4;

%% iterate over differSEIR T1 fits
for ii=1: length(mrQ.inputdata_seir.IT)
if ii>1  
            display(sprintf(['The SPGR - EPI registration is not good enough \nWe will try another pair of maps.']));
end
    % in each iteration, align the SEIR images with a different image as the target image
    % Keeps track of the variables we use.
    % For details, see inside the function.
    
    Dir_SEIR_Align2Im = fullfile(mrQ.outDir, ['SEIR_Align2Im_' num2str(ii)]);
    SEIR_fit_dir = fullfile(Dir_SEIR_Align2Im, 'fitT1_GS');
    
    
% %     T1ref=mrQ.SEIRfits{ii}.SEIR_epi_T1file;
     T1mov=SPGR_T1_file;
     refBMFile=mrQ.SEIRfits{ii}.SEIR_epi_Maskfile;
   % maybe using the skull stripped map is better
     refBM=readFileNifti(refBMFile);refBM=logical(refBM.data);
     T1seirFile=mrQ.SEIRfits{ii}.SEIR_epi_T1file;
     T1seir=readFileNifti(T1seirFile);
     T1seir.data(~refBM)=0;
     xform=T1seir.qto_xyz;  
     [T1path,T1name]=fileparts(T1seirFile);
     T1ref=[T1path,'/skulstripped_',T1name,'.gz'];
    dtiWriteNiftiWrapper(T1seir.data,xform,T1ref);
    %% create a temp folder for the ...
    
    FSL_Path = fullfile(SEIR_fit_dir,'FSL');
    if ~exist(FSL_Path,'dir'), mkdir(FSL_Path); end
    
    %% register the images
    %% linear regitration to being the image close together:
    
    affMat=fullfile(FSL_Path,'linReg.mat');
    [~,T1movName]=fileparts(SPGR_T1_file);
    alignedT1Mov=fullfile(FSL_Path,['Aligned_',T1movName,'.gz']);
    flrt_cmd=['flirt -ref ',T1ref,' -in ',T1mov,' -dof 12 -omat ',affMat,' -out ',alignedT1Mov];
    system(flrt_cmd);
  
        
    %% nonlinear registration
    nonLinMat=fullfile(FSL_Path,'non_linReg');
     % the commented command uses the aligned image when registering spgr to
    % epi. thiscreate double blurring of the image. instead, we try toi
    % only use the linear registration as an initial guess for the
    % nonlinear registration. 
    %   fnrt_cmd=['fnirt  --ref=',T1ref,' --in=',alignedT1Mov,' --cout=',nonLinMat];
    fnrt_cmd=['fnirt  --ref=',T1ref,' --in=',alignedT1Mov,' --aff=',affMat,' --cout=',nonLinMat];
    system(fnrt_cmd);


    %% apply warp
    
    T1_spgr_epi=fullfile(FSL_Path,['Warped_',T1movName,'.gz']);
    % the commented command uses the aligned image when registering spgr to
    % epi. thiscreate double blurring of the image. instead, we try toi
    % only use the linear registration as an initial guess for the
    % nonlinear registration. 
    % warp_cmd=['applywarp --ref=',T1ref,' --in=',alignedT1Mov,' --warp=',nonLinMat,' --out=',T1_spgr_epi];
    
        warp_cmd=['applywarp --ref=',T1ref,' --in=',alignedT1Mov,' --warp=',nonLinMat,' --premat=',affMat,' --out=',T1_spgr_epi];

    system(warp_cmd);
    
    
    
    %% check the registration
    
    MovingScaleConstat=1000;% to compare SPGR T1 in sec to SEIR T1 in ms.
    [mrQ.Ants_Info.QuantAntsScore]=mrQ_QuantAnts(mrQ.SEIRfits{ii}.SEIR_epi_T1file,T1_spgr_epi,MovingScaleConstat);
    
    %% We will keep the files only if they are the current best Ants registration
    
%     if this is the best registration yet...
    if min(mrQ.Ants_Info.QuantAnts2High(:)) > mrQ.Ants_Info.QuantAntsScore
        
        
        mrQ.Ants_Info.WARP_SPGR_EPI = nonLinMat;
        mrQ.Ants_Info.T1_spgr_epi=T1_spgr_epi;
        
        % Saving the epi alignment parameters
        mrQ.Ants_Info.SEIR_SPGR_Curent_AlignNums=[ii,jj];
        mrQ. Ants_Info.SEIR_SPGR_Curent_AlignDirs={SEIR_fit_dir, spgr_initDir};
        mrQ. Ants_Info.SPGRFieldName=SPGRFieldName;
        
        %saving the SEIR files
        mrQ.SEIR_epi_T1file = mrQ.SEIRfits{ii}.SEIR_epi_T1file;
        mrQ.SEIR_epi_resnormfile = mrQ.SEIRfits{ii}.SEIR_epi_resnormfile;
        mrQ.SEIR_epi_fitFile = mrQ.SEIRfits{ii}.SEIR_epi_fitFile;
        mrQ.SEIR_epi_M0file = mrQ.SEIRfits{ii}.SEIR_epi_M0file;
        mrQ.SEIR_epi_Maskfile =mrQ.SEIRfits{ii}.SEIR_epi_Maskfile;
    end
    
    %         check if the current (best) registration is good enough
    if mrQ.Ants_Info.QuantAntsScore< mrQ.QuantAntsThresh;
        mrQ.Ants_Info.QuantAnts2High(ii,jj)=0;
        save(mrQ.name,'mrQ');
        
        break
    else
        % if the registration quality did not pass the threshold, we will
        % keep the values of the registration quality and move on to the
        % next iteration.
        
        mrQ.Ants_Info.QuantAnts2High(ii,jj)=mrQ.Ants_Info.QuantAntsScore;
        save(mrQ.name,'mrQ');
        
    end
end
