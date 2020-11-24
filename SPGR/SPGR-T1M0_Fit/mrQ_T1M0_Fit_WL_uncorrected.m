function mrQ=mrQ_T1M0_Fit_WL_uncorrected(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
%
% In this function, the T1-M0 WL fit is performed for the uncorrected B1. 
%
% ~INPUTS~
%        mrQ: The mrQ structure
%     B1File: Location of the uncorrected B1 file (NIfTI)
%   MaskFile: Location of the mask file (NIfTI)
%     outDir: Directory to save the data
%    dataDir: Directory for the spgr data
%    clobber: Overwrite existing data and reprocess. [Default = false]
%
% ~OUTPUTS~
%        mrQ: The updated mrQ structure
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%


%%  LSQ or LINEAR fit of M0 and T1:
%  Use the SunGrid to accelerate this fit

%% I. Check INPUTS and set defaults

if notDefined('dataDir');
    dataDir = mrQ.InitSPGR.spgr_initDir;
end

if notDefined('outDir');
    outDir =dataDir;
end

% Clobber flag: 
% Overwrite existing fit, if it already exists, and redo the T1 fit
if notDefined('clobber')
    clobber = false;
end

B1=niftiRead(B1File);
B1=double(B1.data);



if isfield(mrQ,'FullMaskFile')
    MaskFile=mrQ.FullMaskFile;
end

if notDefined('MaskFile')
    MaskFile= fullfile(outDir,'FullMask.nii.gz');
end



    
LWfit = 1;


%% II. Load aligned data

outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);

%% V. Weighted linear fit
if LWfit==1

    T1WLFfile= fullfile(outDir,['T1_map_Wlin_unCorrected.nii.gz']);
    M0WLFfile= fullfile(outDir,['M0_map_Wlin_unCorrected.nii.gz']);
    
    if (exist( T1WLFfile,'file') && exist( M0WLFfile,'file')  && ~clobber),
        
        disp(['Loading existing T1 and M0 weighted linear fit'])
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
        
        
        [T1L,M0L] = relaxFitT1(cat(4,s(:).imData),flipAngles,tr,B1);

        % Compute a linear fit of the the T1 estimate for all voxels.
        % M0: PD = M0 * G * exp(-TE / T2*).
        if  exist(MaskFile,'file')
            HeadMask=readFileNifti(MaskFile);
            HeadMask=logical(HeadMask.data);
        else
            [HeadMask, mrQ]=CalculateFullMask(mrQ,T1L,M0L,outDir); %see below for function
        end
        Gain=double(HeadMask);
        
%         [T1w, T1,M0w, MO] = mrQ_T1M0_LWFit(s,HeadMask,tr,flipAngles,Gain,B1,outDir,xform,SunGrid,[],sub,mrQ.proclus);
        
%         [T1w, T1,M0w, M0] = mrQ_T1M0_LWFit(s,HeadMask,tr,flipAngles,Gain,B1,outDir,xform,mrQ.SunGrid);
                [T1w, T1,M0w, M0] = mrQ_T1M0_LWFit(s,HeadMask,tr,flipAngles,Gain,B1,outDir,xform,mrQ);

        % Save the T1 and PD data
        dtiWriteNiftiWrapper(single(T1w), xform,T1WLFfile);
        dtiWriteNiftiWrapper(single(M0w), xform, M0WLFfile);
        
        mrQ.T1_LWFit_unCorrected=T1WLFfile;
        mrQ.M0_B1_LWFit_unCorrected=M0WLFfile;
        
    end
end

       save(mrQ.name,'mrQ');

       %
      %%%
     %%%%%
    %%%%%%%
     %%%%%
      %%%
       %
end      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
function [mask, mrQ]=CalculateFullMask(mrQ,t1,M0,outDir)
%% one more mask anywhere we have signal

if isfield(mrQ,'LinFit')
    HM_path=mrQ.LinFit.HeadMask;
    HM = readFileNifti(mrQ.LinFit.HeadMask); xform= HM.qto_xyz;    HM=logical(HM.data);
    BM=  readFileNifti(mrQ.LinFit.BrainMask); BM=logical(BM.data);
else
    HM_path=mrQ.HeadMask;
    HM = readFileNifti(mrQ.HeadMask);     xform= HM.qto_xyz;  HM=logical(HM.data);
    BM=  readFileNifti(mrQ.BrainMask); BM=logical(BM.data);
end

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


% Create a white-matter mask using a range of T1 values within the original
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
end

% make a head mask that definitely includes the brain mask
mask=logical(mask1+mask +HM+BM);

FullMaskFile= fullfile(outDir,'FullMask.nii.gz');
dtiWriteNiftiWrapper(single(mask), xform,FullMaskFile);

% If SNR is bad use the head mask instead of the full mask:
if isfield(mrQ,'SNR')
    if strcmp(mrQ.SNR,'low') 
        FullMaskFile=HM_path;
        mask=HM;
        warning(['Low SNR - using the head mask instead of the full mask']);
    end
end    
mrQ.FullMaskFile=FullMaskFile;
end
