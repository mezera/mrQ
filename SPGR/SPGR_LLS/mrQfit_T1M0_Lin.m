function [mrQ]=mrQfit_T1M0_Lin(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
%  [mrQ]=mrQfit_T1M0_Lin(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
%
% INPUT:
%       mrQ:        the mrQ structure. 
%       B1File:     nifti file of a B1 map, default is ones
%                   matrice
%       MaskFile:   nifti file of brain mask if not defined it will be
%                   created using mrAnatExtractBrain, and folowed by some
%                   corections. 

%       clobber: 
% 
% OUTPUT:
%                   the function creates nifti filesin the outDir, and saves their
%                   location in the mrQ sructure which is the output of the
%                   function. this function will create niftis for the (if
%                   don't exist): headmask, brainmask, and linearly fitted
%                   T1.
% 
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
%
% EXAMPLE USAGE:
%    *** add relevant example ***
%    [mrQ.AnalysisInfo]=mrQfit_T1M0_ver2(mrQ.spgr_initDir,mrQ.lsq,mrQ.SEIRepiDir,mrQ.SPGR_coilWeights,mrQ.runfreesurfer,mrQ.sub);
%
% (c) Stanford University, VISTA Lab

% Author: Aviv Mezer  date 01.18.2011
% Copyright, Stanford University 2011
% rewrite by AM. June 16 2011
% rewrite by AM. July 29 2011 %fit for the B1 fitting to be based on local
% regression


%%

%% I. Check INPUTS and set defaults



if notDefined('dataDir');
    dataDir = mrQ.spgr_initDir;
end
if notDefined('outDir');
    outDir =dataDir;
end

% Clobber flag. Overwrite existing fit if it already exists and redo the T1
%  fit
if notDefined('clobber')
    clobber = false;
end


%% II. Load aligned data


outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);




%% IV. Setup and save AnalysisInfo file.
% % a structure that keep records of the T1 fit
% infofile=fullfile(outDir,'AnalysisInfo.mat');
% if (exist(infofile,'file'))
%     load (infofile);
% end
%
%
% AnalysisInfo.fieldStrength      = fieldStrength;
% AnalysisInfo.outDir             = outDir;
% AnalysisInfo.dataDir            = dataDir;
% AnalysisInfo.complexFlag        = complexFlag;
% AnalysisInfo.lsqfit             = lsqfit;
% AnalysisInfo.LWfit             = LWfit;
% AnalysisInfo.SEIRepi_Dir        = SEIRepi_Dir;
% AnalysisInfo.T1MOfitdata        = date;
% AnalysisInfo.sub                = sub;
% save(infofile,'AnalysisInfo');
%


%% V. Calculate initial fits without correction

% !!!!!!!!!!!!!!!!!!!!!!!
% All the inital parts can be done in lower resolution. Should save it and
% use it for the B1 gain part.

if notDefined('B1File')
    
    disp('Initial fits with no correction are now being calculated...');
    B1 = ones(size(s(1).imData));
else
    B1=niftiRead(B1File);
    B1=double(B1.data);
end



%% VI. Linear fit to estimate T1 and M0 no B1 (this will be used to fit B1)



t1file   = fullfile(outDir,['T1_LFit.nii.gz']);
M0file   = fullfile(outDir,['M0_LFit.nii.gz']);

% Read in existing T1, M0 and the brain mask (if they exist)
if~( exist(t1file,'file') &&  exist(M0file,'file')  && ~clobber),
   % disp([' Loding existing  T1  M0 and brain mask' ]);
        
    
    % Now we fit T1 and M0:
    disp('1. Performing linear fit of T1 and M0');
    
    % Specify the flip angle and TR: s2 is loaded when the dat_aligned.mat
    % file is loaded above.
    flipAngles = [s(:).flipAngle];
    tr         = [s(:).TR];
    
    % Check that all TRs are the same.
    if ~all(tr == tr(1))
        error('TR''s do not match!');
    end
    tr = tr(1);
    
    %
    [t1,M0]     = relaxFitT1(cat(4,s(:).imData),flipAngles,tr,B1);
    
    t1_copy=t1;M0_copy=M0; % make a copy
    
    if notDefined('MaskFile')
        HMfile   = fullfile(outDir,'HeadMask.nii.gz');  % brain mask
        BMfile   = fullfile(outDir,'brainMask.nii.gz'); %head mask
        % Create the brain mask
        [brainMask,checkSlices] = mrAnatExtractBrain(M0, mmPerVox, 0.5,outDir);
      %  eval(['! rm ' outDir '/bet* '])
        
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
        
%         SAVE head mask and brain mask as niftis
        dtiWriteNiftiWrapper(single(HM), xform, HMfile);
        dtiWriteNiftiWrapper(single(brainMask), xform, BMfile);
           mrQ.HeadMask=HMfile;
    mrQ.BrainMask=BMfile;
    else
        
        brainMask=niftiRead(MaskFile);
        brainMask=logical(brainMask.data);
        
    end
    
     
    t1(~brainMask) = 0;
    M0(~brainMask) = 0;
    
    % Save the T1 and PD data
    dtiWriteNiftiWrapper(single(t1), xform, t1file);
    dtiWriteNiftiWrapper(single(M0), xform, M0file);
    
    mrQ.T1_LFit=t1file;
    mrQ.M0_LFit_HM=M0file;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if ~notDefined('HMfile')
        t1=t1_copy; M0=M0_copy;
        t1(~HM) = 0; 
        M0(~HM) = 0;
        
        t1fileHM = fullfile(outDir,['T1_LFit_HM.nii.gz']);
        M0fileHM = fullfile(outDir,'M0_LFit_HM.nii.gz');
        
        
        dtiWriteNiftiWrapper(single(t1), xform, t1fileHM);
        dtiWriteNiftiWrapper(single(M0), xform, M0fileHM);
        mrQ.T1_LFit_HM=t1fileHM;
        mrQ.M0_LFit_HM=M0fileHM;
    
    end
    
end;




    
