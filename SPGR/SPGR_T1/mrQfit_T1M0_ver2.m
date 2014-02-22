function [AnalysisInfo]=mrQfit_T1M0_ver2(mrQ,clobber)
%
%   mrQfit_T1M0_ver2(dataDir,lsqfit,SEIRepi_Dir,outDir,complexFlag,sub)
%
% SPGR data with a range of flip angles are re-sampled and the HLF is
% estimated.
%
% INPUTS:
%       dataDir            - The directory where the aligned SPGR exists
%                            named dat_aligned.mat
%
%       lsqfit             - [1] - for lsqfit of T1 using the Grid
%                            (recommended).
%                            [0] - for linear fit (a fast default)
%
%       SEIRepi_Dir:       - Path to the epi SEIR directory.
%                            The epi SEIR directory should contain the
%                            out-put from fit mrQ_fitSEIR_T1 (see above)
%       coilWeights         - if a coill Weighting data exsist we will use it . in tje other case we
%                              will use the coil Weighting as given from the magnets
%
%      runfreesurfer        - if a freesurfer is needed we will make a call
%                           for that after we crate a sysntetic T1-Weighted image
%
%       clobber:              Overwrite existing data and reprocess. [default = false]
%
%       outDir             - Path to where you want the data to be
%                            written to. If the out directory is same then
%                            the data directory default, leave it empty.
%
%
%
%       complexFlag:       - if SEIRepi ws complex data set to 1.
%                            default = zero
%
%       sub                - Subject name used when SGE cmd is executed.
%
%
% OUTPUT:
%       AnalysisInfo       -The function  tries to keep track of the last
%                           runs so it will save a structure name AnalysisInfo  in the outDir.
%                           AnalysisInfo.mat keeps track of the input and dates of the
%                           file in the directory. That file will be updated when this function
%                           is run again and when the other relevant functions are used like
%                           mrQ_PD and mrQ_HLF.
%                           maps like T1 M0 B1, etc. are written to the
%                           output directory ('outDir').
%
%
%
%
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
%
% EXAMPLE USAGE:
%    [mrQ.AnalysisInfo]=mrQfit_T1M0_ver2(mrQ.spgr_initDir,mrQ.lsq,mrQ.SEIRepiDir,mrQ.SPGR_coilWeights,mrQ.runfreesurfer,mrQ.sub);
%
% (c) Stanford University, VISTA Lab

% Author: Aviv Mezer  date 01.18.2011
% Copyright, Stanford University 2011
% rewrite by AM. June 16 2011
% rewrite by AM. July 29 2011 %fit for the B1 fitting to be based on local
% regression

%#ok<*FNDSB>
%#ok<*ASGLU>
%#ok<*NASGU>
%%

%% I. Check INPUTS and set defaults

if ~isfield(mrQ,'spgr_initDir');
    mrQ.spgr_initDir = pwd;
end
dataDir = mrQ.spgr_initDir;


if ~isfield(mrQ,'lsq');
    mrQ.lsq = 1;
end
lsqfit = mrQ.lsq;

if ~isfield(mrQ,'SEIRepiDir')
    error('can not fit SEIR epi to spgr the SEIRepi_Dir path is missing')
end


if ~isfield(mrQ,'outDir');
    mrQ.outDir = mrQ.spgr_initDir;
end
outDir = mrQ.outDir;


if ~isfield(mrQ,'complexFlag');
    mrQ.complexFlag = 0;
end
complexFlag = mrQ.complexFlag;

if ~isfield(mrQ,'runfreesurfer')
    mrQ.runfreesurfer=0;
end
runfreesurfer = mrQ.runfreesurfer;

if ~isfield(mrQ,'sub');
    % This is a job name we get from the for SGE
    [~, mrQ.sub] = fileparts(fileparts(fileparts(fileparts(fileparts(dataDir)))));
    disp([' Subject name for lsq fit is ' mrQ.sub]);
end
sub = mrQ.sub;
if ~isfield(mrQ,'SPGR_coilWeight');
    mrQ.SPGR_coilWeight = 0;
end
coilWeights = mrQ.SPGR_coilWeight;

if ~isfield(mrQ,'name');
    mrQ.name = fullfile(outDir,'mrQParams');
end

if ~isfield(mrQ,'lsq');
    mrQ.lsq = 1;
end
lsqfit = mrQ.lsq;

if ~isfield(mrQ,'SunGrid');
    mrQ.SunGrid = 1;
end
SunGrid = mrQ.SunGrid;

save(mrQ.name,'mrQ')

% Clobber flag. Overwrite existing fit if it already exists and redo the T1
% B1 fit

if notDefined('clobber')
    clobber = false;
end





%% I. Check INPUTS and set defaults

if (~exist('dataDir','var') || isempty(dataDir)),
    dataDir = pwd;
end

if (~exist('lsqfit','var') || isempty(lsqfit)),
    lsqfit = 0;
end

if (~exist(mrQ.SEIRepiDir,'dir') || isempty(mrQ.SEIRepiDir)),
    error('can not fit SEIR epi to spgr the SEIRepi_Dir path is missing')
else
    SEIRepiDir = mrQ.SEIRepiDir;
    SEIRepi_Dir = mrQ.SEIRepiDir;
end

if (~exist('outDir','var') || isempty(outDir)),
    outDir = dataDir;
end

if(~exist(outDir,'dir')), mkdir(outDir); end

if (~exist('complexFlag','var') || isempty(complexFlag)),
    complexFlag = 0;
end

if notDefined('sub')
    % This is a job name we get from the for SGE
    [~, sub] = fileparts(fileparts(fileparts(fileparts(fileparts(dataDir)))));
    disp([' Subject name for lsq fit is ' sub]);
end

if notDefined('coilWeights')
    coilWeights = 1;
end

if notDefined('runfreesurfer')
    runfreesurfer = 0;
end

% Clobber flag. Overwrite existing fit if it already exists and redo the T1
% B1 fit
if notDefined('clobber')
    clobber = false;
end


%% II. Load aligned data

if coilWeights==1
    outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data
    outFile1 = fullfile(dataDir,'dat_alignedBest.mat');%with coilWeights data
else
    outFile  = fullfile(dataDir,'dat_aligned.mat'); %in this case the is no  coilWeights data
    outFile1 = fullfile(dataDir,'dat_aligned.mat');
end


disp(['Loading aligned data from ' outFile '...']);

load(outFile);
load(outFile1);

if notDefined('s2')
    s2=s;
end
%% III. Get field strength
if isfield(s,'fieldStrength');

    fieldStrength = s(1).fieldStrength;
    disp(['Magnet field strength is: ' num2str(fieldStrength) 'T']);
else
    fieldStrength = input(['please provid the Magnet Filed  (0.5 1.5 or 3)  ']) ;
end

if  (fieldStrength~=(1.5) && fieldStrength~=(3) && fieldStrength~=(0.5))
    error('please provid a Magnet Filed  (0.5 1.5 or 3)')
end


%% IV. Setup and save AnalysisInfo file.
% a structure that keep records of the T1 fit
infofile=fullfile(outDir,'AnalysisInfo.mat');
if (exist(infofile,'file'))
    load (infofile);
end


AnalysisInfo.fieldStrength      = fieldStrength;
AnalysisInfo.outDir             = outDir;
AnalysisInfo.dataDir            = dataDir;
AnalysisInfo.complexFlag        = complexFlag;
AnalysisInfo.lsqfit             = lsqfit;
AnalysisInfo.SEIRepi_Dir        = SEIRepi_Dir;
AnalysisInfo.T1MOfitdata        = date;
AnalysisInfo.sub                = sub;
save(infofile,'AnalysisInfo');



%% V. Calculate initial fits without correction

% !!!!!!!!!!!!!!!!!!!!!!!
% All the inital parts can be done in lower resolution. Should save it and
% use it for the B1 gain part.
disp('Initial fits with no correction are now being calculated...');
B1 = ones(size(s(1).imData));



%% VI. Linear fit to estimate T1 and M0 no B1 (this will be used to fit B1)

t1fileHM = fullfile(outDir,['T1_LFit_HM.nii.gz']);
HMfile   = fullfile(outDir,'HeadMask.nii.gz');
M0fileHM = fullfile(outDir,'M0_LFit_HM.nii.gz');

BMfile   = fullfile(outDir,'brainMask.nii.gz');
t1file   = fullfile(outDir,['T1_LFitB.nii.gz']);
M0file   = fullfile(outDir,['M0_LFitB.nii.gz']);

% Read in existing T1, M0 and the brain mask (if they exist)
if( exist(t1file,'file') &&  exist(M0file,'file') &&  exist(BMfile,'file') && ~clobber),
    disp([' Loding existing  T1  M0 and brain mask' ]);

    t1=readFileNifti(t1file);
    t1=t1.data;

    M0=readFileNifti(M0file);
    M0=M0.data;

    brainMask = readFileNifti(BMfile);
    brainMask=logical(brainMask.data);
else


    % Now we fit T1 and M0:
    disp('1. Performing linear fit of T1 and M0');

    % Specify the flip angle and TR: s2 is loaded when the dat_aligned.mat
    % file is loaded above.
    flipAngles = [s2(:).flipAngle];
    tr         = [s2(:).TR];

    % Check that all TRs are the same.
    if ~all(tr == tr(1))
        error('TR''s do not match!');
    end
    tr = tr(1);

    % Fitting routine from RFD and NS: s2 is the aligned data from
    % mrQ_initSPGR.m and loaded above in cell II.
    % Using the data from the aligned and aligned best
    % files here and fitting twice. (s and s2).
    [t1,M0]     = relaxFitT1(cat(4,s2(:).imData),flipAngles,tr,B1);
    flipAngles0 = [s(:).flipAngle];
    [tt,M01]    = relaxFitT1(cat(4,s(:).imData),flipAngles0,tr,B1);

    % Create the brain mask
    [brainMask,checkSlices] = mrAnatExtractBrain(M01, mmPerVox, 0.5,outDir);
    eval(['! rm ' outDir '/bet* '])
    % Replace all nan values in the brain mask with zeros.
    for dd=1:length(s)
        brainMask(isnan(s(dd).imData))=0;
    end;
    %findind the first and the last slice of the brain
    [~, ~, z]=ind2sub(size(brainMask),find(brainMask));

    % Create the head mask the head mask make the registration with the
    % SEIR imags better using ANts softwere

    % we find the voxel that have signal geater then 2std below the
    % mean brain signal
    cutV=mean(s(1).imData(brainMask)) -2*std(s(1).imData(brainMask));
    %noise estraction - we will need to genralize this somehow espesialy

    % selecting the relevant slices (this is not beutiful)
    HM=B1;
    HM(s(1).imData<cutV)=0;
    HM(:,:,1:min(z)+3)=0;
    HM(:,:,max(z)+3:end)=0;
    for dd=1:length(s)
        HM(isnan(s(dd).imData))=0;
    end;

    % feeling the holes in the mask.
    for i=1:size(HM,3)
        HM(:,:,i)=imfill(HM(:,:,i),'holes');
    end;
    HM=logical(HM);

    % Save the masks,T1 and M0 (PD) data
    dtiWriteNiftiWrapper(single(HM), xform, HMfile);
    dtiWriteNiftiWrapper(single(brainMask), xform, BMfile);

    tt(~HM) = 0;
    M01(~HM) = 0;
    t1(~brainMask) = 0;
    M0(~brainMask) = 0;
    %t1(t1>5) = 5;

    % Save the T1 and PD data
    dtiWriteNiftiWrapper(single(t1), xform, t1file);
    dtiWriteNiftiWrapper(single(M0), xform, M0file);
    dtiWriteNiftiWrapper(single(tt), xform, t1fileHM);
    dtiWriteNiftiWrapper(single(M01), xform, M0fileHM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end;

%% VII. FIT B1
%  Required Inputs: SEIRepi_T1  AlignFile
% we will Align the SPGR and the SEIR by linear Warp using ANTS softwere

B1file=fullfile(outDir,['B1_Map.nii.gz']);
AlignFile=fullfile(outDir,['SEIRepiSPGRAlign_best_RB.mat']);

% If B1 has already been computed we load it, if not then we fit it.
if exist(B1file,'file') && ~clobber
    B1=readFileNifti(B1file);
    B1=double(B1.data);

else
    %%

    if complexFlag ==0
        SET1file=[SEIRepi_Dir '/fitT1_GS/T1FitNLSPR_SEIR_Dat_T1.nii.gz'];
    elseif complexFlag ==1
        SET1file=[SEIRepi_Dir '/fitT1_GS/T1FitNLS_SEIR_Dat_T1.nii.gz'];

    end;
    if exist(SET1file,'file')
    else
        [file_se path_se]= uigetfile(pwd,'Select the T1 SEIR map');
        SET1file=[path_se file_se];
    end



    SET1Fitfile=[SEIRepi_Dir '/fitT1_GS/T1FitNLSPR_SEIR_Dat.mat'];
    SET1=readFileNifti(SET1file);
    SE_Xform=SET1.qto_xyz;
    pixdim=SET1.pixdim;
    clear SET1


    B1epifile=fullfile(outDir,['B1epi_RB_full_best.nii.gz']);
    B1epiResidfile=fullfile(outDir,['ResidB1epi_NLW_full_best.nii.gz']);

    % If the fit in SEIR space has already been computed we load it
    if exist(B1epifile,'file') && ~clobber

    else

        % Load the SEIR_to SPGR Aligned "best" file - or create it
        if exist(AlignFile,'file') && ~clobber
            load (AlignFile)
        else
            flipAngles = [s2(:).flipAngle];

            if ~isfield(mrQ,'AligndSPGR'); % if the Align image was not crated we will make them know
                for j=1:length(s2)
                    ref   = fullfile(outDir,['Align' num2str(flipAngles(j)) 'deg']);
                    dtiWriteNiftiWrapper(single(s2(j).imData), xform,ref);
                    mrQ.AligndSPGR{j}=ref;

                end
                save(mrQ.name,'mrQ');
            end

            [AnalysisInfo,Res]=mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile,mrQ.AligndSPGR);

        end


        %%% FIT B1 by lsq fit compare T1 SEIR(Res{1}.im) to the multi flip
        % angle Res{3:end}.im % USE sge make the fit faster

        intM0= double(mean(M0(brainMask)));
        flipAngles = [s2(:).flipAngle];
        tr = [s2(:).TR];

        if(~all(tr == tr(1))), error('TR''s do not match!'); end
        tr = tr(1);

        % Load or Create the tissuemask from all non-zero points
        tisuuemaskFile=fullfile(outDir,['maskepiF.nii.gz']);

        if exist(tisuuemaskFile,'file') && ~clobber
            tisuuemask_=readFileNifti(tisuuemaskFile);
            tisuuemask_=logical(tisuuemask_.data);

        else

            tisuuemask_ =zeros(size(Res{1}.im));

            % Binarize the mask
            tisuuemask_(find(Res{1}.im>10 & Res{3}.im>0))=1;

            % Create a logical array from the tissue mask to index the
            % non-zero locations
            tisuuemask_=logical(tisuuemask_);

            % Save the tissue mask
            dtiWriteNiftiWrapper(single(tisuuemask_), SE_Xform, tisuuemaskFile);

        end;

        if clobber && (exist([outDir '/tmpSG'],'dir'))
            % in the case we start over and there are old  fits we will
            % deleat them
            eval(['! rm -r ' outDir '/tmpSG']);
        end

        % USE sge make the B1 fit faster


        % This is lsq fit that uses the grid but you can make it not use
        % SGE: see help inside mrQ_fitB1_LSQ
        [B1 resNorm dd] = mrQ_fitB1_LSQ(Res, tisuuemask_, tr,flipAngles, outDir, intM0, SE_Xform, SunGrid, 1, [sub 'B1fit'],mrQ.proclus);
        dtiWriteNiftiWrapper(single(B1), SE_Xform, B1epifile);
        dtiWriteNiftiWrapper(single(resNorm), SE_Xform, B1epiResidfile);

    end

    %%% Move B1 from epi space to SPGR smooth and upsample

    % keyboard

    % Smooth the SEIR B1 map by local and then fill the gap by global
    % regretions and then register back the smoothed B1 to the SPGR space using ANTS
    [B1,AnalysisInfo]=mrQ_smmothL_B1(B1epifile,AlignFile,outDir,B1file,xform,[],AnalysisInfo,t1fileHM,SET1Fitfile,B1epiResidfile);


    % Create the synthetic T1 weighted images and save them to disk
    AnalysisInfo.T1wSynthesis=mrQ_T1wSynthesis(dataDir,B1file,outDir);
    save(infofile,'AnalysisInfo');



    if runfreesurfer==1
        subjID=['freesufer_' AnalysisInfo.sub];
        [AnalysisInfo.freeSufer_subdir]=mrQ_Callfs_autosegment(subjID, AnalysisInfo.T1wSynthesis);
        save(infofile,'AnalysisInfo');
    end


end


%% VIII. LSQ or LINEAR fit of M0 and T1:
%  Use the sun-grid to excelerate this fit

% make a head mask that include for sure the brain mask
HM = readFileNifti(HMfile);

HeadMask=logical(HM.data +brainMask);

%%% LSQ FIT
if lsqfit==1,
    %% lsq fit of M0 and T1 use the sun-grid to excelerate this fit
    disp('Fitting T1 and PD by lsq: This takes time - SGE can be used!!!');
    T1lsqfile= fullfile(outDir,['T1_map_lsq.nii.gz']);
    M0lsqfile= fullfile(outDir,['M0_map_lsq.nii.gz']);

    % If the fits have been performed then we simply load them.

    if (exist( T1lsqfile,'file') && exist( M0lsqfile,'file') && ~clobber),

        disp(['loading exsisting T1 and M0 lsq fit'])
        T1=readFileNifti(T1lsqfile);
        M0=readFileNifti(M0lsqfile);
        T1=double(T1.data);
        M0=double(M0.data);

    else


        if clobber && (exist([outDir '/tmpSG'],'dir'))
            % in the case we start over and there are  old fits we will
            % deleat them
            eval(['! rm -r ' outDir '/tmpSG']);
        end



        disp(['Fiting lsq T1 and M0']);
        flipAngles = [s2(:).flipAngle];
        tr = [s(:).TR];

        Gain=double(HeadMask);

        % LSQ fit of M0 and T1: Use the sun-grid to excelerate this fit
        [T1,M0] = mrQ_fitT1PD_LSQ(s2,HeadMask,tr,flipAngles,M0,t1,Gain,B1,outDir,xform,SunGrid,[],sub,mrQ.proclus);

        %

        % Save the T1 and M0 data
        dtiWriteNiftiWrapper(single(T1), xform,T1lsqfile);
        dtiWriteNiftiWrapper(single(M0), xform, M0lsqfile);

        AnalysisInfo.T1lsqfile=T1lsqfile;

    end

    %%% LINEAR FITTING (lsqfit ~=1) Linear fit is used to calculate T1 and M0
    %%% (Linear fitting can bias the fit but it's very fast)

end
disp([' linear fits of T1 and PD !!!'] );
T1LFfile= fullfile(outDir,['T1_map_lin.nii.gz']);
M0LFfile= fullfile(outDir,['M0_map_lin.nii.gz']);
%

if (exist( T1LFfile,'file') && exist( M0LFfile,'file')  && ~clobber),

    disp(['loading exsisting T1 and M0 linear fit'])
    T1=readFileNifti(T1LFfile);
    M0=readFileNifti(M0LFfile);
    T1=double(T1.data);
    M0=double(M0.data);

else

    disp('Performing linear fit of T1 and M0...');
    flipAngles = [s2(:).flipAngle];
    tr = [s(:).TR];

    % Check that all TRs are the same across all the scans in S
    if(~all(tr == tr(1))), error('TR''s do not match!'); end
    tr = tr(1);

    % Compute a linear fit of the the T1 estimate for all voxels.
    % M0: PD = M0 * G * exp(-TE / T2*).
    [T1,M0] = relaxFitT1(cat(4,s(:).imData),flipAngles,tr,B1);

    % Zero-out the values that fall outside of the brain mask
    T1(~HeadMask) = 0;
    M0(~HeadMask) = 0;


    % Save the T1 and PD data
    dtiWriteNiftiWrapper(single(T1), xform,T1LFfile);
    dtiWriteNiftiWrapper(single(M0), xform, M0LFfile);

    AnalysisInfo.T1LFfile=T1LFfile;

end;

save(infofile,'AnalysisInfo');






