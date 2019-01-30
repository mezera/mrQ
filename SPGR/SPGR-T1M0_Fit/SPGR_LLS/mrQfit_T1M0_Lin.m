function [Outputs]=mrQfit_T1M0_Lin(mrQ,dataDir,outDir,B1File,MaskFile,clobber)


%  [mrQ]=mrQfit_T1M0_Lin(mrQ,B1File,MaskFile,outDir,dataDir,clobber)
%
% This function creates NIfTI files in the outDir, and saves their
% location in the updated mrQ structure, which is the output of the
% function. This function will create NIfTIs for the following (if they
% don't exist): headmask, brainmask, and linearly fitted T1.
%
% INPUT:
%       mrQ:        The mrQ structure.
%       B1File:     NIfTI file of a B1 map. Default is a ones matrix.
%       MaskFile:   NIfTI file of brain mask. If not defined, it will be
%                   created using mrAnatExtractBrain.m, and followed by some
%                   corrections.
%       outDir:     Directory to which the data will be saved.
%       dataDir:    Directory of the data.
%       clobber:    Overwrite existing data and reprocess. [Default = false]
%
% OUTPUT:
%       mrQ:        The mrQ structure, updated.
%
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
%
% EXAMPLE USAGE:
%    *** add relevant example ***
%    [mrQ.AnalysisInfo]=mrQfit_T1M0_ver2(mrQ.spgr_initDir,mrQ.lsq,mrQ.SEIRepiDir,mrQ.SPGR_coilWeights,mrQ.runfreesurfer,mrQ.sub);
%
% (C) Stanford University, VISTA Lab
%
% Author: Aviv Mezer, 01.18.2011
% rewrite by AM. June 16 2011
% rewrite by AM. July 29 2011
% rewrite by AM. May 29 2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016

%     %fit for the B1 fitting to be based on local regression


%% I. Check INPUTS and set defaults

if notDefined('dataDir');
    dataDir = mrQ.mrQ.InitSPGR.spgr_initDir ;
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


%% III. Setup and save AnalysisInfo file.
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


%% IV. Calculate initial fits without correction

% !!!!!!!!!!!!!!!!!!!!!!!
% All the initial parts can be done in lower resolution. Should save it and
% use it for the B1 gain part.

if notDefined('B1File')
    
    disp('Initial fits with no correction are now being calculated...');
    B1 = ones(size(s(1).imData));
else
    B1=niftiRead(B1File);
    B1=double(B1.data);
end



%% V. Linear fit to estimate T1 and M0 no B1 (this will be used to fit B1)

t1file   = fullfile(outDir,['T1_LFit.nii.gz']);
M0file   = fullfile(outDir,['M0_LFit.nii.gz']);

% Read in existing T1, M0 and the brain mask (if they exist)
%if~( exist(t1file,'file') &&  exist(M0file,'file')  && ~clobber),
    % disp([' Loding existing  T1  M0 and brain mask' ]);
    
    
    % Now we fit T1 and M0:
    disp('1. Performing linear fit of T1 and M0');
    
    % Specify the flip angle and TR:
    % s2 is loaded when the dat_aligned.mat file is loaded above.
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
    
   % if notDefined('MaskFile')
        HMfile   = fullfile(outDir,'HeadMask.nii.gz');  % brain mask
        BMfile   = fullfile(outDir,'brainMask.nii.gz'); %head mask
        % Create the brain mask
       
        [brainMask,checkSlices] = mrAnatExtractBrain(M0, mmPerVox, 0.5,outDir);
         out = fullfile(tempdir,'bet_tmp');
        delete([out,'.img'])
        delete([out,'.hdr'])
        %  eval(['! rm ' outDir '/bet* '])
        
        % Replace all nan values in the brain mask with zeros.
        for dd=1:length(s)
            brainMask(isnan(s(dd).imData))=0;
        end;
        
        %find the first and the last slice of the brain
        [~, ~, z]=ind2sub(size(brainMask),find(brainMask));
        
        %% VI. Create the head mask.
        % The head mask makes the registration with the SEIR images better
        % using ANTs software
        
        % Find the voxels that are more than 2 stdev below the
        % mean brain signal. First we will chose the image with max signal.
        % this approch not working allways so we add some checks.
       
        Zrange=median(z)-3:median(z)+3;
        NBM=length(find(brainMask(:,:,Zrange))); %number of voxel in the center of the brain mask
        NB1=length(find(B1(:,:,Zrange))); %number of voxel in the image in the same slices
        for ii=1:length(s)
            M(ii)= mean(s(ii).imData(brainMask));
        end
        
        
        
        jj=find(M==max(M)); % max mean signal.
        cutV=mean(s(jj).imData(brainMask)) -2*std(s(jj).imData(brainMask));
        
        if cutV<0 || cutV<nanmedian(s(jj).imData(~brainMask)) 
            % if our cut value (cutV) for no head voxels is negative or bellow the
            % median non brain values it is probably wrong.
                warning('   Selecting a threshold for head-mask is challenging, we will try a different strategy. Please,be aware, inspection the head-mask-output is advised.');

            cutV=nanmedian(s(jj).imData(~brainMask))*2; % we can try an alternative: the double of the non brain values
%            An alternative solution is to use a different jj value. We
%            didn't implament that. it can be done later or manually.
        end
        %noise extraction - we will need to generalize this somehow especially
        
        % selecting the relevant slices (this is not beautiful)
        HM=B1;
        %%
        HM(s(jj).imData<cutV)=0;
        %AM SF test the effect of comment the below two lines

       % HM(:,:,1:min(z)+3)=0;
       % HM(:,:,max(z)+3:end)=0;
        %%
        for dd=1:length(s)
            HM(isnan(s(dd).imData))=0;
        end
        
        % filling the holes in the mask
        for i=1:size(HM,3)
            HM(:,:,i)=imfill(HM(:,:,i),'holes');
        end;
        
        NHM=length(find(HM(:,:,Zrange))); %number of voxel in the center of the head mask

        if NHM/NB1>0.9 || NBM/NHM<0.6 || NBM/NHM>1 % let's check that the head mask is reasnable. biger than the brain mask ans smaller than the image size.
            warning(' The head-mask seem wrong. Please check! Please consider to use the brain-mask as the head-mask or edit the head-mask manually. This may introduce error in PD fitting.')
        end
        
        %%
        HMtmp=HM;

                HMtmp(t1<0.4)=0;

        for i=1:size(HM,3)
           FracNonBrainT1(i)= length(find(HMtmp(:,:,i)))./ length(find(HM(:,:,i))) ;
        end;
        %%
        HM=logical(HM);
        
        %%  VII. SAVE head mask and brain mask as niftis
        dtiWriteNiftiWrapper(single(HM), xform, HMfile);
        dtiWriteNiftiWrapper(single(brainMask), xform, BMfile);
        mrQ.HeadMask=HMfile;
        mrQ.BrainMask=BMfile;
   % else
        
%         brainMask=niftiRead(MaskFile);
%         brainMask=logical(brainMask.data);
%         
%     end
    
    t1(~brainMask) = 0;
    M0(~brainMask) = 0;
    
    %% VIII. SAVE the T1 and PD data
    dtiWriteNiftiWrapper(single(t1), xform, t1file);
    dtiWriteNiftiWrapper(single(M0), xform, M0file);
    
%     mrQ.T1_LFit=t1file;
%     mrQ.M0_LFit=M0file;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if ~notDefined('HMfile')
        t1=t1_copy; M0=M0_copy;
        t1(~HM) = 0;
        M0(~HM) = 0;
        
        t1fileHM = fullfile(outDir,['T1_LFit_HM.nii.gz']);
        M0fileHM = fullfile(outDir,'M0_LFit_HM.nii.gz');
        
        
        dtiWriteNiftiWrapper(single(t1), xform, t1fileHM);
        dtiWriteNiftiWrapper(single(M0), xform, M0fileHM);
%         mrQ.T1_LFit_HM=t1fileHM;
%         mrQ.M0_LFit_HM=M0fileHM;
%         mrQ.FracNonBrainT1=FracNonBrainT1;

   % end
    % mrQ.HeadMask=HMfile;
    %         mrQ.BrainMask=BMfile;
    %
    %     mrQ.T1_LFit=t1file;
    %     mrQ.M0_LFit=M0file;
    %  mrQ.T1_LFit_HM=t1fileHM;
    %         mrQ.M0_LFit_HM=M0fileHM;
    %         mrQ.FracNonBrainT1=FracNonBrainT1;
    
%end;

Outputs.HeadMask=HMfile;
Outputs.BrainMask=BMfile;
Outputs.T1_LFit=t1file;
Outputs.M0_LFit=M0file;
Outputs.T1_LFit_HM=t1fileHM;
Outputs.M0_LFit_HM=M0fileHM;
Outputs.FracNonBrainT1=FracNonBrainT1;
