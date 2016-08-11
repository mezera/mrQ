 function [opt]=mrQ_resamp4G_fit(opt,outMm,interp,T1reg)
% [opt]=mrQ_resamp4G_fit(opt,outMm,interp,T1reg)
%
% This function resamples the M0 (4D), T1 (3D) and brain mask (3D) images
% to the outMm resolution and saves those file locations to the opt
% structure so they can be used for Gain Fit. 
%
%
%   ~INPUTS~
%           opt: Structure with all the parameters.
%         outMm: The resample (undersample) resolution of those images.
%                   (Default is 2mm x 2mm x 2mm.)
%        interp: Defines the interpolation method according to the SPM
%                    standard, as used in mrAnatResliceSpm.m
%                    Default is 1.
%         T1reg:
%
%  ~OUTPUTS~
%          opt: The updated opt structure.
%
%

% use mrQ_UpSampG4PD_fit.m
%M0file='/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1/AligncombineCoilsM0.nii.gz'

%% I. Set default values
 if notDefined('outMm')
     outMm=[2 2 2];
 end
 
 if notDefined('interp')
     interp = 1;
 end
 
 opt.Resamp=1;
 
%% II. Multi-coil M0

M0=readFileNifti(opt.M0file);

sz=size(M0.data);
bb = mrAnatXformCoords(M0.qto_xyz, [1 1 1; sz(1:3)]);

if length(sz)>3
    for ii=1:sz(4)
        [M0UnderSamp(:,:,:,ii), M0UnderSamp_Xform] = mrAnatResliceSpm(double(M0.data(:,:,:,ii)), inv(M0.qto_xyz), bb, outMm, interp, 0);
        
    end
    sos=sqrt(sum(M0UnderSamp.^2,4));
    clear M0
else
    [M0UnderSamp, M0UnderSamp_Xform] = mrAnatResliceSpm(double(M0.data), inv(M0.qto_xyz), bb, outMm, interp, 0);
    sos=M0UnderSamp;
end
    
filename=fullfile(opt.outDir,['M0resmp_' num2str(outMm(1)) '_' num2str(outMm(2)) '_' num2str(outMm(3)) '.nii.gz']);
dtiWriteNiftiWrapper(single(M0UnderSamp),M0UnderSamp_Xform,filename);

opt.M0file_Org=opt.M0file;
opt.M0file=filename;
 clear  M0UnderSamp M0UnderSamp_Xform filename

%% III. Brain mask
BM=readFileNifti(opt.BMfile);
sz=size(BM.data);

%bb = mrAnatXformCoords(BM.qto_xyz, [1 1 1; sz(1:3)]);

[BMUnderSamp, BMUnderSamp_Xform] = mrAnatResliceSpm(double(BM.data(:,:,:)), inv(BM.qto_xyz), bb, outMm, interp, 0);

%% IV. Tissue mask
if isfield(opt,'TissueFile')
    TM=readFileNifti(opt.TissueFile);
    sz=size(BM.data);
    %bb = mrAnatXformCoords(BM.qto_xyz, [1 1 1; sz(1:3)]);
    [TissueMUnderSamp, TissueMUnderSamp_Xform] = mrAnatResliceSpm(double(TM.data(:,:,:)), inv(TM.qto_xyz), bb, outMm, interp, 0);
end

%% V. T1
if isfield(opt,'T1file')
    T1=readFileNifti(opt.T1file);
    
    
    [T1UnderSamp, T1UnderSamp_Xform] = mrAnatResliceSpm(double(T1.data(:,:,:)), inv(T1.qto_xyz), bb, outMm, interp, 0);
    
    
    
% Clean and save T1 and BM
    % Clean
    
    % These values are probably out of the original brian mask and are a
    % result of the resampling process
    
    T1UnderSamp((T1UnderSamp<0.4))=0; % brain tisuue
    T1UnderSamp((T1UnderSamp>5))=5; % too high T1
    BMUnderSamp(T1UnderSamp==0)=0; % not brain by T1 values
    BMUnderSamp(T1UnderSamp==5)=0; % not brain by T1 values
    
    if isfield(opt,'TissueFile')
        TissueMUnderSamp(T1UnderSamp==0)=0; % not brain by T1 values
        TissueMUnderSamp(T1UnderSamp==5)=0; % not brain by T1 values
    end
    
end

BMUnderSamp(BMUnderSamp<1)=0; %not sure this is brain any more

B=ones(3,3,3); % a 3X3x3 mask

%to be sure, let's only keep the voxels that are brain voxels
C=convn(BMUnderSamp,B,'same'); 

BMUnderSamp(C<9)=0;
BMUnderSamp(sos<=0)=0;
BMUnderSamp=logical(BMUnderSamp);

if isfield(opt,'TissueFile')
    TissueMUnderSamp(C<9)=0;
    TissueMUnderSamp(sos<=0)=0;
    TissueMUnderSamp=logical(TissueMUnderSamp);
end

    % Save
%BM
filename=fullfile(opt.outDir,['BMresmp_' num2str(outMm(1)) '_' num2str(outMm(2)) '_' num2str(outMm(3)) '.nii.gz']);
dtiWriteNiftiWrapper(single(BMUnderSamp),BMUnderSamp_Xform,filename);
opt.BMfile_Org=opt.BMfile;
opt.BMfile=filename;

if isfield(opt,'TissueFile')
    filename=fullfile(opt.outDir,['TMresmp_' num2str(outMm(1)) '_' num2str(outMm(2)) '_' num2str(outMm(3)) '.nii.gz']);
    dtiWriteNiftiWrapper(single(TissueMUnderSamp),BMUnderSamp_Xform,filename);
    opt.TissueFile_Org=opt.TissueFile;
    opt.TissueFile=filename;
end

%clear BM BMUnderSamp BMUnderSamp_Xform filename

if isfield(opt,'T1file')
    %T1
    filename=fullfile(opt.outDir,['T1resmp_' num2str(outMm(1)) '_' num2str(outMm(2)) '_' num2str(outMm(3)) '.nii.gz']);
    dtiWriteNiftiWrapper(single(T1UnderSamp),T1UnderSamp_Xform,filename);
    
    opt.T1file_Org=opt.T1file;
    opt.T1file=filename;
end

%clear T1 T1UnderSamp T1UnderSamp_Xform filename

