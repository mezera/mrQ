function mrQ_WF(outDir,PDfile,CSFfile)
%
% mrQ_WF(outDir,PDfile,CSFfile)
%
% Compute the water fraction: Divide PS by the CSF to get WF
%
% This function will take a freesurfer segmentation in nifti format
% (aparc+aseg.nii.gz) and the directory of the data analysis and a
% M0file and T1 images.
%
% It will generate a smooth white matter mask for for the coil gain
% estimatin wiht polynoyal with a degree (defult 3) . and estimate it and
% save a clean PD and WF maps if M0 is not provided it will look for
% PD_lsqnabs_SEIR_cfm.nii.gz in the outDir
%
% If you like to see the freesurfer on T1 defult yes
%
% INPUT:
%   outDir  -  Directory containing the aligned SPGR data.
%   PDfile  -  PDsqrt_fitboxMedian_1.nii.gz
%   CSFfile -  csf_FS_T1.nii.gz
%
% OUTPUT:
%    WFfile = fullfile(outDir,'WF_fitG2n.nii.gz');
%    WFfile = fullfile(outDir,'WF_fit.nii.gz');
% 
%
% (C) Stanford University, VISTA
%


%% Check INPUTS


if(exist('CSFfile','var') &&  ~isempty(CSFfile) )
    disp(['Loading CSF data from ' CSFfile '...']);
else
    CSFfile = fullfile(outDir,'csf_FS_T1.nii.gz');    
    if(exist(CSFfile,'file'))
        disp(['Loading CSF data from ' CSFfile '...']);
    else
        CSFfile = mrvSelectFile('r','*.nii.gz','Select CSF File',outDir);
        if isempty(CSFfile)
            error('User Cancelled.');
        end
    end
end


if(exist('PDfile','var') &&  ~isempty(PDfile))
    disp(['Loading PD data from ' PDfile '...']);
else
    PDfile = fullfile(outDir,'PDsqrt_fitboxMedian_1.nii.gz');    
    if(exist(PDfile,'file'))
        disp(['Loading PD data from ' PDfile '...']);
    else
        PDfile = mrvSelectFile('r','*.nii.gz','Select PD File',outDir);
        if isempty(PDfile)
            error('User Cancelled.');
        end
    end
end


%% Load CSF & PD data, the Analysis Info and Brain Mask


CSF = readFileNifti(CSFfile);
CSF = double(CSF.data);

PD = readFileNifti(PDfile);
PD = double(PD.data);


infofile = fullfile(outDir,'AnalysisInfo.mat');
load(infofile);

AnalysisInfo.CSFfile  = CSFfile;
AnalysisInfo.PDfile   = PDfile;
%AnalysisInfo.M0forWF = M0file;
AnalysisInfo.CSFdate  = date;
save(infofile,'AnalysisInfo');

BMfile = fullfile(outDir,'brainMask.nii.gz');

if ~exist(BMfile,'file')
    BMfile = mrvSelectFile('r','*.nii.gz','Select Brain Mask File',outDir);
    if isempty(BMfile)
        error('User Cancelled.');
    end
end

disp(['Loading brain Mask data from ' BMfile '...']);
brainMask = readFileNifti(BMfile);
xform     = brainMask.qto_xyz;
mmPerVox  =  brainMask.pixdim;
brainMask = logical(brainMask.data);


%%

% The PD was fit as the sqrt of the data when the gain was fit
PD = PD.^2; 

% The CSF is not in the WM values and not double it (may be wrong for kids)
CSF1 = CSF & PD > 7000 & PD < 14000;

% Kernel smoothing density estimate
[d dd] = ksdensity( PD(find(CSF1)), [min(PD(find(CSF1))):0.01:max(PD(find(CSF1)))] );

%  CalibrationVal = mean(PD(find(CSF)));
CalibrationVal = dd(find(d==max(d)));

% median(PD(find(CSF)));

% Compute WF
WF = PD./CalibrationVal;

% Threshold WF between 0 and 1.
WF(WF>1) = 1;  WF(WF<0) = 0;

% Save WF 
WFfile = fullfile(outDir,'WF_fit.nii.gz');
dtiWriteNiftiWrapper(single(WF), xform, WFfile);


%% Load tissue file and T1 fit

% Load the freesurfer tissue file (where is this generated???)
FSfile = fullfile(outDir,'FS_tissue.nii.gz');
if ~exist(FSfile,'file')
    FSfile = mrvSelectFile('r','*.nii.gz','Select Freesurfer tissue file',outDir);
    if isempty(FSfile)
        error('User cancelled.');
    end
end

fs     = readFileNifti(FSfile);
fs     = fs.data;

% T1 fit. This is loaded, but does not appear to be used. 
T1file = fullfile(outDir,'T1_lsq_GLr.nii.gz');
if ~exist(T1file,'file')
    T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',outDir);
    if isempty(T1file)
        warning('No T1 fit file was selected.')
    end
end

T1 = readFileNifti(T1file);
T1 = T1.data;


%%

% KNK code: Use polynomial basis functions to fit PD.  we try to fit the
% elements marked with a 1 in (fs==2), then 2 in (fs==2), and so forth. We
% enforce the constraint that different cases have the same surface shape
% but allow each case to have its own gain factor.  see
% constructpolynomialmatrix3d.m for details on the basis functions.
[params,gains,rs] = fit3dpolynomialmodel(PD,(fs==2),2);
Imsz1 = size(PD);

[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),2);

G = reshape(Poly1*params(:),Imsz1);

PD1 = zeros(size(brainMask));
PD1(brainMask) = PD(brainMask)./G(brainMask);


% The CSF is not in the WM values and not double it (may be wrong for kids)
CSF1 = CSF & PD1>1 & PD1< 2; 


% Kernel smoothing density estimate
[d dd] = ksdensity(PD1(find(CSF1)), [min(PD1(find(CSF1))):0.001:max(PD1(find(CSF1)))] );

CalibrationVal = dd(find(d==max(d))); % median(PD(find(CSF))); %  CalibrationVal = mean(PD(find(CSF)));

WF1 = PD1./CalibrationVal;

% Threshold the WF
WF1(WF1>1) = 1;  WF1(WF1<0)=0;

% Write the WF file to disk
WFfile = fullfile(outDir,'WF_fitG2n.nii.gz');
dtiWriteNiftiWrapper(single(WF1), xform, WFfile);

return



 
%% Old code

%     G2file=fullfile(outDir,['G2.nii.gz']);
%     dtiWriteNiftiWrapper(single(G), xform,G2file);

%PolFitonData(fs==2,PD,100)
%     [params,gains,rs] = fit3dpolynomialmodel(PD,(fs==2),4);
%     Imsz1=size(PD);
%     [Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),4);
%
%     G1 = reshape(Poly1*params(:),Imsz1);
%
%              G4file=fullfile(outDir,['G4.nii.gz']);
%     dtiWriteNiftiWrapper(single(G1), xform,G4file);
%
