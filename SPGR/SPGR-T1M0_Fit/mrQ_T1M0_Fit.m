function mrQ=mrQ_T1M0_Fit(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
% mrQ=mrQ_T1M0_Fit(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
%
% Different types of fit of T1 to M0: linear, least squares (LSQ), linear weighted
%
% ~INPUTS~
%        mrQ: The mrQ structure
%     B1File: Location of the B1 file (NIfTI)
%   MaskFile: Location of the mask file (NIfTI)
%     outDir: Directory to save the data
%    dataDir: Directory for the spgr data
%    clobber: Default is false 
%
% ~OUTPUTS~
%        mrQ: The updated mrQ structure
%
%


%%  LSQ or LINEAR fit of M0 and T1:
%  Use the SunGrid to accelerate this fit

%% I. Check INPUTS and set defaults

if notDefined('dataDir');
    dataDir = mrQ.spgr_initDir;
end

if notDefined('outDir');
    outDir =dataDir;
end

% Clobber flag: 
% Overwrite existing fit, if it already exists, and redo the T1 fit
if notDefined('clobber')
    clobber = false;
end

if isfield(mrQ,'B1FileName')
    B1File=mrQ.B1FileName;
end

if notDefined('B1File')
    disp('Initial fits with no correction are now being calculated...');
    B1 = ones(size(s(1).imData));
else
    B1=niftiRead(B1File);
    B1=double(B1.data);
end

if isfield(mrQ,'FullMaskFile')
    MaskFile=mrQ.FullMaskFile;
end

if notDefined('MaskFile')
       MaskFile= fullfile(outDir,'FullMask.nii.gz');
end

if ~isfield(mrQ,'lsq');
    mrQ.lsq = 1;
end

lsqfit = mrQ.lsq;

if ~isfield(mrQ,'LW');
    mrQ.LW = false;
end

if lsqfit==0 
    mrQ.LW = true;
end

LWfit = mrQ.LW;


%% II. Load aligned data

outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);

%% III. Linear Fit 
% Linear fit is used to calculate T1 and M0.
% (Linear fitting can bias the fit, but it's very fast)

disp([' Linear fits of T1 and PD !!!']);
T1LFfile= fullfile(outDir,['T1_map_lin.nii.gz']);
M0LFfile= fullfile(outDir,['M0_map_lin.nii.gz']);

if (exist( T1LFfile,'file') && exist( M0LFfile,'file')  && exist( MaskFile,'file') && ~clobber),
    
    disp(['Loading existing T1 and M0 linear fit'])
    
    T1L=readFileNifti(T1LFfile);
    M0L=readFileNifti(M0LFfile);
    HeadMask=readFileNifti(MaskFile);
    
    T1L=double(T1L.data);
    M0L=double(M0L.data);
    HeadMask=logical(HeadMask.data);
        
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
    
    % let's calculate a new mask 
    [HeadMask, mrQ]=CalculateFullMask(mrQ,T1L,M0L,outDir); %see below for function
    % Zero-out the values that fall outside of the brain mask
    T1L(~HeadMask) = 0;
    M0L(~HeadMask) = 0;    
    
    % Save the T1 and PD data
    dtiWriteNiftiWrapper(single(T1L), xform,T1LFfile);
    dtiWriteNiftiWrapper(single(M0L), xform, M0LFfile);
    
    mrQ.T1_B1_LFit=T1LFfile;
    mrQ.M0_B1_LFit=M0LFfile;
end;


if ~isfield(mrQ,'FullMaskFile')
        [HeadMask, mrQ]=CalculateFullMask(mrQ,T1L,M0L,outDir);
end


%% IV. LSQ Fit
if lsqfit==1,
    % LSQ fit of M0 and T1. Use the SunGrid to accelerate this fit
    
    disp('Fitting T1 and PD by LSQ: This takes time - SGE can be used!!!');
    T1lsqfile= fullfile(outDir,['T1_map_lsq.nii.gz']);
    M0lsqfile= fullfile(outDir,['M0_map_lsq.nii.gz']);
    
    % If the fits have been performed, then we simply load them.
    
    if (exist( T1lsqfile,'file') && exist( M0lsqfile,'file') && ~clobber),
        
        disp(['Loading existing T1 and M0 LSQ fit'])
        T1=readFileNifti(T1lsqfile);
        M0=readFileNifti(M0lsqfile);
        T1=double(T1.data);
        M0=double(M0.data);
        
    else
        if clobber && (exist([outDir '/tmpSG'],'dir'))
            % in the case we start over and there are old fits, 
            % we will delete them
            eval(['! rm -r ' outDir '/tmpSG']);
        end      
        
        disp(['Fitting LSQ T1 and M0']);
        flipAngles = [s2(:).flipAngle];
        tr = [s(:).TR];
        
        Gain=double(HeadMask);
        
        % LSQ fit of M0 and T1: Use the SunGrid to accelerate this fit
        [T1,M0] = mrQ_fitT1PD_LSQ(s2,HeadMask,tr,flipAngles,M0L,T1L,Gain,B1,outDir,xform,SunGrid,[],sub,mrQ.proclus);
              
        % Save the T1 and M0 data
        dtiWriteNiftiWrapper(single(T1), xform,T1lsqfile);
        dtiWriteNiftiWrapper(single(M0), xform, M0lsqfile);
        
        mrQ.T1_B1_lsqFit=T1lsqfile;
        mrQ.M0_B1_lsqFit=M0lsqfile;
    end
end

%% V. Weighted linear fit
if LWfit==1
    disp([' Linear fits of T1 and PD !!!'] );
    T1WLFfile= fullfile(outDir,['T1_map_Wlin.nii.gz']);
    T1LFfile= fullfile(outDir,['T1_map_lin.nii.gz']);
    M0LFfile= fullfile(outDir,['M0_map_Wlin.nii.gz']);
    M0WLFfile= fullfile(outDir,['M0_map_Wlin.nii.gz']);
    
    if (exist( T1WLFfile,'file') && exist( M0WLFfile,'file')  && ~clobber),
        
        disp(['Loading existing T1 and M0 linear fit'])
        T1=readFileNifti(T1WLFfile);
        M0=readFileNifti(M0WLFfile);
        T1=double(T1.data);
        M0=double(M0.data);
        
    else
        
        disp('Performing weighted linear fit of T1 and M0 with B1 correction...');
        flipAngles = [s(:).flipAngle];
        tr = [s(:).TR];
        
        % Check that all TRs are the same across all the scans in S
        if(~all(tr == tr(1))), error('TR''s do not match!'); end
        tr = tr(1);
        
        % Compute a linear fit of the the T1 estimate for all voxels.
        % M0: PD = M0 * G * exp(-TE / T2*).
        
        Gain=double(HeadMask);
        
%         [T1w, T1,M0w, MO] = mrQ_T1M0_LWFit(s,HeadMask,tr,flipAngles,Gain,B1,outDir,xform,SunGrid,[],sub,mrQ.proclus);
        
        [T1w, T1,M0w, M0] = mrQ_T1M0_LWFit(s,HeadMask,tr,flipAngles,Gain,B1,outDir,xform,mrQ.SunGrid);

        % Save the T1 and PD data
        dtiWriteNiftiWrapper(single(T1), xform,T1LFfile);
        dtiWriteNiftiWrapper(single(M0), xform, M0LFfile);
        dtiWriteNiftiWrapper(single(T1w), xform,T1WLFfile);
        dtiWriteNiftiWrapper(single(M0w), xform, M0WLFfile);
        
        mrQ.T1_B1_LWFit=T1WLFfile;
        mrQ.M0_B1_LWFit=M0WLFfile;
        
        mrQ.T1_B1_LFit=T1LFfile;
        mrQ.M0_B1_LFit=M0LFfile;
        
    end
end

       save(mrQ.name,'mrQ');

function [mask, mrQ]=CalculateFullMask(mrQ,t1,M0,outDir)
%% one more mask anywhere we have signal

HM = readFileNifti(mrQ.HeadMask);     HM=logical(HM.data);
BM=  readFileNifti(mrQ.BrainMask); BM=logical(BM.data);

% M=mean(t1(BM));
% S=std(t1(BM));
% 
% up=min(10000,M+2*S);
% down=max(0,M-2*S);
% t1(t1<down)=down;
% t1(t1>up)=up;
% 
% t1(isnan(t1))=down;
% 
% t1(isinf(t1))=up;
% 
% mask1=t1>prctile(t1(BM),0.1) & t1<up;
% [mask1] = ordfilt3D(mask1,6);
% for i=1:size(mask1,3)
%     mask1(:,:,i)=imfill(mask1(:,:,i),'holes');
% end;
% %%%%
% Create a white-matter mask using a range of t1 values within the original
% brain mask

wmmask    = BM & t1>0.850 & t1<1.050;

% Perform 3-D order-statistic filtering on M0
[wmmask1] = ordfilt3D(wmmask,6);
wmmask    = wmmask & wmmask1; clear wmmask1;

[params1,gains,rs] = fit3dpolynomialmodel(M0,wmmask==1,2);
wmmask             = logical(wmmask);

 Imsz = size(BM);

% Construct a three-dimenstional polynomial matrix - KNK code
[Poly,str] = constructpolynomialmatrix3d(Imsz,find(ones(Imsz)),2);
%   Returns <Poly>, a matrix of dimensions length(find(ones(Imsz))) x N
%   with polynomial basis functions evaluated at <find(ones(Imsz))> in
%   the columns.  the polynomial basis functions are evaluated
%   over the range [-1,1] which is presumed to correspond to
%   the beginning and ending element along each of the three dimensions.
%   (if a dimension has only one element, the values are all set to 1.)
%   also, return <str>, the algebraic expression that corresponds to
%   the columns of <Poly>.  'x' refers to the first matrix dimension; 'y'
%   refers to the second matrix dimension; 'z' refers to the third
%   matrix dimension.

% Calculate the gain and proton density
Gain = reshape(Poly*params1',Imsz);
M0   = M0./Gain;


% mask

M=mean(M0(BM));
S=std(M0(BM));

up=min(10000,M+3*S);
down=max(0,M-3*S);
M0(M0<down)=down;
M0(M0>up)=up;

M0(isnan(M0))=down;

M0(isinf(M0))=up;

mask=M0>prctile(M0(BM),0.1) & M0<up;
[mask] = ordfilt3D(mask,6);
for i=1:size(mask,3)
    mask1(:,:,i)=imfill(mask(:,:,i),'holes');
end;

mask=logical(mask1+mask +HM+BM);


% make a head mask that definitely includes the brain mask
FullMaskFile= fullfile(outDir,'FullMask.nii.gz');
dtiWriteNiftiWrapper(single(mask), mrQ.xform,FullMaskFile);
mrQ.FullMaskFile=FullMaskFile;

