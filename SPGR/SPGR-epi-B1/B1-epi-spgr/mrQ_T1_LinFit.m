function mrQ=mrQ_T1_LinFit(mrQ,B1File,MaskFile,outDir,dataDir)

%% I. Check INPUTS and set defaults

if notDefined('dataDir');
    dataDir = mrQ.Ants_Info.SEIR_SPGR_Curent_AlignDirs{2};
end

if notDefined('outDir');
    outDir =dataDir;
end

% Clobber flag:
% Overwrite existing fit, if it already exists, and redo the T1 fit

if notDefined('B1File')   
    if isfield(mrQ,'B1FileName')
        B1File=mrQ.B1FileName;
    end   
end
if isfield(mrQ,'HeadMask')
    MaskFile=mrQ.HeadMask;
end


%% II. Load aligned data

outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);

% % % % % % % % % % %
if notDefined('B1File')
    disp('Initial fits with no correction are now being calculated...');
    B1 = ones(size(s(1).imData));
else
    B1=niftiRead(B1File);
    B1=double(B1.data);
end

% % % % % % % % % % %
if notDefined('MaskFile')
    HeadMask= ones(size(s(1).imData));
else
    HeadMask=readFileNifti(MaskFile);
    HeadMask=logical(HeadMask.data);
end

%% III. Linear Fit
% Linear fit is used to calculate T1 and M0.
% (Linear fitting can bias the fit, but it's very fast)

disp([' Linear fits of T1 and PD !!!']);
T1LFfile= fullfile(outDir,['T1_lin_unCorrected.nii.gz']);
M0LFfile= fullfile(outDir,['M0_lin_unCorrected.nii.gz']);

if (exist( T1LFfile,'file') && exist( M0LFfile,'file'))
    
    disp(['Loading existing T1 and M0 linear fit'])
    
    T1L=readFileNifti(T1LFfile);
    M0L=readFileNifti(M0LFfile);
    
    T1L=double(T1L.data);
    M0L=double(M0L.data);
    
else
    
    disp('Performing linear fit of T1 and M0...');
    
    flipAngles = [s(:).flipAngle];
    tr = [s(:).TR];
    
    % Check that all TRs are the same across all the scans in S
    if(~all(tr == tr(1))), error('TR''s do not match!'); end
    tr = tr(1);
    
    % Compute a linear fit of the the T1 estimate for all voxels.
    % M0: PD = M0 * G * exp(-TE / T2*).
    [T1L,M0L] = relaxFitT1(cat(4,s(:).imData),flipAngles,tr,B1);
    
    % Zero-out the values that fall outside of the brain mask
    T1L(~HeadMask) = 0;
    M0L(~HeadMask) = 0;
    
    % Save the T1 and PD data
    dtiWriteNiftiWrapper(single(T1L), xform,T1LFfile);
    dtiWriteNiftiWrapper(single(M0L), xform, M0LFfile);
    
    mrQ.T1_B1_LFit_unCorrected=T1LFfile;
    mrQ.M0_B1_LFit_unCorrected=M0LFfile;
    
    save(mrQ.name,'mrQ');
end;


