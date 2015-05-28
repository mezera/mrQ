function mrQ_fitT1M0(dataDir,lsqfit,SEIRepi_Dir,outDir,inputAlignFile_old,complexFlag,sub)
% 
% mrQ_fitT1M0(dataDir,lsqfit,SEIRepi_Dir,outDir,inputAlignFile_old,complexFlag,sub)
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
% 
%       outDir             - Path to where you want the data to be
%                            written to. If the out directory is same then
%                            the data directory default, leave it empty.
% 
%       inputAlignFile_old - if the SPGR and SEIR_EPI were aligned before
%                            (using our code), you can specify a .mat file
%                            to be used as a starting point for a new
%                            alignment.
% 
%       complexFlag:       - if SEIRepi ws complex data set to 1. 
%                            default = zero
% 
%       sub                - Subject name used when SGE cmd is executed.
% 
% 
% OUTPUT:
%       This is empty, but maps like T1 M0 B1, etc. are written to the
%       output directory ('outDir'). 
%       
%       The function also tries to keep track of the last
%       run so it will save an AnalysisInfo file in the outDir. 
%       AnalysisInfo.mat keeps track of the input and dates of the
%       file in the directory. That file will be updated when this function
%       is run again and when the other relevant functions are used like
%       mrQ_PD and mrQ_HLF.
%
% 
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
% 
%
% EXAMPLE USAGE:
%       mrQ_fitT1M0(dataDir,1,SEIRepi_Dir,[],[])
%       mrQ_fitT1M0('/biac2/wandell2/data/WMDevo/adult/QuantitativeImaging/104_JY/20110610_0536/SPGR_1/Align0.9375_0.9375_1/',1,'/biac2/wandell2/data/WMDevo/adult/QuantitativeImaging/104_JY/20110610_0536/SEIR_epi_1',[],[])
% 
% 
% (c) Stanford University, VISTA Lab

% Author: Aviv Mezer  date 01.18.2011 
% Copyright, Stanford University 2011 
% rewrite by AM. June 16 2011 
% rewrite by AM. July 29 2011 %fit for the B1 fitting to be based on local
% regression

%#ok<*FNDSB>
%#ok<*ASGLU>
%#ok<*NODEF>
%#ok<*NASGU>

%% I. Check INPUTS and set defaults

if (~exist('dataDir','var') || isempty(dataDir)),
    dataDir = pwd;
end

if (~exist('lsqfit','var') || isempty(lsqfit)),
    lsqfit = 0;
end

if (~exist('SEIRepi_Dir','var') || isempty(SEIRepi_Dir)),
    error('can not fit SEIR epi to spgr the SEIRepi_Dir path is missing')
end

if (~exist('outDir','var') || isempty(outDir)),
    outDir = dataDir;
end
if(~exist(outDir,'dir')), mkdir(outDir); end

if (~exist('inputAlignFile_old','var') || isempty(inputAlignFile_old)),
    inputAlignFile_old = [];
end

if (~exist('complexFlag','var') || isempty(complexFlag)),
    complexFlag = 0;
end

if (~exist('sub','var') || isempty(sub))
    % This is a job name we get from the for SGE
    [~, b] = fileparts(fileparts(fileparts(fileparts(fileparts(dataDir)))));
    sub = b;
    disp([' Subject name for lsq fit is ' b]);
end


%% II. Load aligned data

outFile  = fullfile(dataDir,'dat_aligned.mat');
outFile1 = fullfile(dataDir,'dat_alignedBest.mat');

disp(['Loading aligned data from ' outFile '...']);

load(outFile);
load(outFile1);


%% III. Get field strength

fieldStrength = s(1).fieldStrength; 
disp(['Magnet field strength is: ' num2str(fieldStrength) 'T']);

if  (fieldStrength~=(1.5) && fieldStrength~=(3) && fieldStrength~=(0.5))
    error('please provid a Magnet Filed  (0.5 1.5 or 3)')
end


%% IV. Setup and save AnalysisInfo file

AnalysisInfo.fieldStrength      = fieldStrength;
AnalysisInfo.outDir             = outDir;
AnalysisInfo.dataDir            = dataDir;
AnalysisInfo.complexFlag        = complexFlag;
AnalysisInfo.lsqfit             = lsqfit;
AnalysisInfo.SEIRepi_Dir        = SEIRepi_Dir;
AnalysisInfo.inputAlignFile_old = inputAlignFile_old;
AnalysisInfo.T1MOfitdata        = date;
% AnalysisInfo.localregressionB1 = localregressionB1;

infofile = fullfile(outDir,'AnalysisInfo.mat');
save(infofile,'AnalysisInfo');


%% V. Calculate initial fits without correction

% !!!!!!!!!!!!!!!!!!!!!!!
% All the inital parts can be done in lower resulotion. Should save it and
% use it for the B1 gain part.
disp('Initial fits with no correction are now being calculated...');
B1 = ones(size(s(1).imData));


%% VI. Linear fit to estimate T1 and M0 no B1 (this will be used to fit B1)

% If this is not the first time running this funciton then these files will
% exist. If they do exist we will simply load them. If they don't exist
% then we will create the brain mask and fit T1 and M0.
BMfile = fullfile(outDir,'brainMask.nii.gz');
t1file = fullfile(outDir,'T1_LFitB.nii.gz');
M0file = fullfile(outDir,'M0_LFitB.nii.gz');

% Read in existing T1, M0 and the brain mask (if they exist)
if exist(t1file,'file') && exist(M0file,'file') && exist(BMfile,'file') 
    disp('Loading existing T1, M0 and Brain Mask...');

    t1 = readFileNifti(t1file);
    t1 = t1.data;
    
    M0 = readFileNifti(M0file);
    M0 = M0.data;
    
    brainMask = readFileNifti(BMfile);
    brainMask = logical(brainMask.data);
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
    [brainMask,checkSlices] = mrAnatExtractBrain(M01, mmPerVox, 0.5);
    
    % Replace all nan values in the brain mask with zeros.
    for dd = 1:length(s)
        brainMask(isnan(s(dd).imData)) = 0;
    end

    % Save out the brain mask. 
    dtiWriteNiftiWrapper(single(brainMask), xform, BMfile);
    
    % Remove data that is outside of the brain mask from T1 and B0
    t1(~brainMask) = 0;
    M0(~brainMask) = 0;
    % t1(t1>5)      = 5;

    % Save the T1 and M0 (PD) data
    dtiWriteNiftiWrapper(single(t1), xform, t1file);
    dtiWriteNiftiWrapper(single(M0), xform, M0file);
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end


%% VII. FIT B1 
%  Required Inputs: SEIRepi_T1 SEraw1 AlignFile inputAlignFile_old

AlignFile = fullfile(outDir,'SEIRepiSPGRAlign_best.mat');
B1file    = fullfile(outDir,'B1_fit_lregGx3.nii.gz');

% B1file=fullfile(outDir,['B1_fm.nii.gz']);
% B1file=fullfile(outDir,['B1_Gx2.nii.gz']); % local regression smooth B1 fitting

% If B1 has already been computed we load it, if not then we fit it. 
if exist(B1file,'file')
    B1 = readFileNifti(B1file);
    B1 = double(B1.data);

else
    
    if complexFlag ==0
        SET1file = fullfile(SEIRepi_Dir,'fitT1_GS','T1FitNLSPR_SEIR_Dat_T1.nii.gz');
    elseif complexFlag ==1
        SET1file = fullfile(SEIRepi_Dir,'fitT1_GS','T1FitNLS_SEIR_Dat_T1.nii.gz');

    end
    
    SET1     = readFileNifti(SET1file);
    SE_Xform = SET1.qto_xyz;

    % B1epifile = fullfile(outDir,['B1epi_fm.nii.gz']);
    B1epifile = fullfile(outDir,'B1_fit_full_best.nii.gz');
    
    % If the fit has already been computed we load it
    if exist(B1epifile,'file')
        B1 = readFileNifti(B1epifile);
        B1 = double(B1.data);
        load (AlignFile)
    else
        % Load the SEIRepi Aligned "best" file - or create it
        if exist(AlignFile,'file')
            load (AlignFile)
        else
            
            % ALIGN the SPGR and the epi (semi manual knk code) we need a
            % dicom directory to get information from the header (SEraw1)
            % and the T1 fitted (T1epi)
            
            if complexFlag == 0     % For Magnitude data fit of SEIR epi
                SEIRepi_T1=[SEIRepi_Dir '/fitT1_GS/T1FitNLSPR_SEIR_Dat.mat'];
            elseif complexFlag == 1 % For complex data fit of SEIR epi
                SEIRepi_T1=[SEIRepi_Dir '/fitT1_GS/T1FitNLS_SEIR_Dat.mat'];
            end
            
            load(SEIRepi_T1)
            T1epi(:,:,:) = ll_T1(:,:,:,1);
            
            epiDatDir = fullfile(SEIRepi_Dir,'data');
            Ni = genpath(epiDatDir);
            
            if(isunix)
                Ni = explode(':',genpath(epiDatDir));
            else
                Ni = explode(';',genpath(epiDatDir));
            end

            % This seems hacky - do we know what directory we're looking
            % for? ***
            SEraw1 = Ni{2};
            clear ll_T1;
            
            % Align the SPGR and SEIT EPIs (semi manual knk code) good luck!!!
            % What is RES???
            OAlignFile = fullfile(outDir,'SEIRepiSPGRAlign.mat');
            Res = mrQ_alignSPGR_SEIRepi(T1epi,SEraw1,t1,s2,pwd,1,AlignFile,inputAlignFile_old,OAlignFile);
            load(AlignFile)

        end

        %%% FIT B1 by lsq fit compare T1 SEIR(Res{1}.im) to the multi flip
        % angle Res{3:end}.im % USE sge make the fit faster
        
        intM0 = double(mean(M0(brainMask)));
        flipAngles = [s2(:).flipAngle];
        tr = [s2(:).TR];
        
        if(~all(tr == tr(1))), error('TR''s do not match!'); end
        tr = tr(1);
        
        % Load or Create the tissuemask from all non-zero points
        tisuuemaskFile = fullfile(outDir,'maskepiF.nii.gz');
        % tisuuemaskFile=fullfile(outDir,['maskepi.nii.gz']);
        if exist(tisuuemaskFile,'file')
            tisuuemask_ = readFileNifti(tisuuemaskFile);
            tisuuemask_ = logical(tisuuemask_.data);
        else
            tisuuemask_ = zeros(size(Res{1}.im));
            % Binarize the mask
            tisuuemask_(find(Res{1}.im>0 & Res{3}.im>0)) = 1; 
            
            % Create a logical array from the tissue mask to index the
            % non-zero locations
            tisuuemask_ = logical(tisuuemask_);
            
            % Save the tissue mask
            dtiWriteNiftiWrapper(single(tisuuemask_), SE_Xform, tisuuemaskFile);
        end
        
        % USE sge make the B1 fit faster
        SunGrid = 1; 

        % This is lsq fit that uses the grid but you can make it not use
        % SGE: see help inside mrQ_fitB1_LSQ
        [B1 resNorm dd] = mrQ_fitB1_LSQ(Res, tisuuemask_, tr,flipAngles, outDir, intM0, SE_Xform, SunGrid, 1, [sub 'B1fit']); 
        dtiWriteNiftiWrapper(single(B1), SE_Xform, B1epifile);
        
    end
    
    %%% Move B1 from epi space to SPGR smooth and upsample
    % [B1_h]=mrQ_fitL_GregB1(B1,SE_Xform,Res,s2,outDir,B1file,pixdim,s,mmPerVox,ttr,xform,degree,matlabpoolF)
    B1 = mrQ_fitL_GregB1(B1, SE_Xform, Res, s2, outDir, B1file, pixdim, s, mmPerVox, ttr, xform, [], 1);
    
    % Create the synthetic T1 weighted images and save them to disk
    mrQ_T1wSynthesis(dataDir,B1file,outDir);
end

%% Debug

% if calcT == 1%this won't stay this is for my debug. (aviv)
%keyboard
% else
%     return
% end;


%% VIII. LSQ or LINEAR fit of M0 and T1: 
%  Use the sun-grid to excelerate this fit

% Give the fit a name
name ='GLr'; 

%%% LSQ FIT
if lsqfit == 1
    disp('Fitting T1 and PD by lsq: This takes time - SGE can be used!!!');
    
    T1lsqfile = fullfile(outDir,['T1_lsq_' name '.nii.gz']);
    M0lsqfile = fullfile(outDir,['M0_lsq_' name '.nii.gz']);
        
    % T1lsqfile= fullfile(outDir,['T1_lsqnabs_' num2str(Slabmm) '.nii.gz']);
    % M0lsqfile= fullfile(outDir,['PD_lsqnabs_' num2str(Slabmm) '.nii.gz']);


    % If the fits have been performed then we simply load them. 
    if (exist(T1lsqfile,'file') && exist(M0lsqfile,'file')),
        
        disp('Loding existing T1 and M0 lsq fit');
        T1 = readFileNifti(T1lsqfile);
        M0 = readFileNifti(M0lsqfile);
        T1 = double(T1.data);
        M0 = double(M0.data);
        % keyboard
        
    else

        disp('Fiting lsq T1 and M0...');
        flipAngles = [s2(:).flipAngle];
        tr = [s(:).TR];

        Gain = double(brainMask);
        
        % LSQ fit of M0 and T1: Use the sun-grid to excelerate this fit
        SunGrid = 1; 
        [T1,M0] = mrQ_fitT1PD_LSQ(s2,brainMask,tr,flipAngles,M0,t1,Gain,B1,outDir,xform,SunGrid,[],sub);

        % Zero-out the values that fall outside of the brain mask
        T1(~brainMask) = 0;
        M0(~brainMask) = 0;

        % Save the T1 and M0 data
        dtiWriteNiftiWrapper(single(T1), xform,T1lsqfile);
        dtiWriteNiftiWrapper(single(M0), xform, M0lsqfile);

    end

    %%%%%%%%%%%

    % %    calculation of PD WF HLF fh SIC from the fitted data
    %     WFlsqfile = fullfile(outDir,['Wf_lsq_' name '.nii.gz']); 
    %          Wf=zeros(size(brainMask));
    %          Wf(brainMask)=M0(brainMask)./Gain(brainMask);
    %          dtiWriteNiftiWrapper(single(Wf), xform, WFlsqfile);
    %     disp(['calculate the lsq HLF  WF R1b and SIC'])
    %
    %
    % %     [HLF,Wf,fh,dd] = hlfCompute(T1,pdC,xform,outDir,MagnetFile,[],0,0);
    %     dtiWriteNiftiWrapper(single(HLF), xform, fullfile(outDir,['HLF_lsq_test1.nii.gz']));
    %     dtiWriteNiftiWrapper(single(Wf), xform, fullfile(outDir,['Wf_lsq_test1.nii.gz']));
    %     dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,['fh_lsq_test1.nii.gz']));
    %
    %     [Rb, SIC] = AperntBondT1(T1,brainMask,Wf,MagnetFile);
    %     dtiWriteNiftiWrapper(single(SIC), xform, fullfile(outDir,['SIC_lsq_test1.nii.gz']));
    %         dtiWriteNiftiWrapper(single(Rb), xform, fullfile(outDir,['R1b_lsq_test1.nii.gz']));
    
    %%%%%%%%%
    
    
%%% LINEAR FITTING (lsqfit ~=1) Linear fit is used to calculate T1 and M0
%%% (Linear fitting can bias the fit but it's very fast)
else 
    disp('Performing linear fits of T1 and PD !!!');
    
    T1LFfile= fullfile(outDir,['T1_lin' name '.nii.gz']);
    M0LFfile= fullfile(outDir,['PD_lin_' name '.nii.gz']);
    
    % T1LFfile= fullfile(outDir,['T1_LF_B1_l.nii.gz']);
    % M0LFfile= fullfile(outDir,['PD_LF_B1_l.nii.gz']);

    if (exist( T1LFfile,'file') && exist( M0LFfile,'file') ),

        disp('loding exsist T1 and M0 linear fit')
        T1 = readFileNifti(T1LFfile);
        M0 = readFileNifti(M0LFfile);
        T1 = double(T1.data);
        M0 = double(M0.data);

    else

        disp('Performing linear fit of T1 and M0...');
        flipAngles = [s2(:).flipAngle];
        tr = [s(:).TR];
        
        % Check that all TRs are the same across all the scans in S
        if(~all(tr == tr(1))), error('TR''s do not match!'); end
        tr = tr(1);

        % Compute a linear fit of the the T1 estimate for all voxels.
        % M0: PD = M0 * G * exp(-TE / T2*).
        [T1,M0] = relaxFitT1(cat(4,s2(:).imData),flipAngles,tr,B1);

        % Zero-out the values that fall outside of the brain mask
        T1(~brainMask) = 0;
        M0(~brainMask) = 0;
        
        % Set an upper-thresh on t1
        % t1(t1>5)      = 5;

        % Save the T1 and PD data
        dtiWriteNiftiWrapper(single(T1), xform, T1LFfile);
        dtiWriteNiftiWrapper(single(M0), xform, M0LFfile);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  Calculation of PD WF HLF fh R1b and SIC from the fitted data
    %     pdC=zeros(size(brainMask));
    %     pdC(brainMask)=M0(brainMask)./Gain(brainMask);
    %     dtiWriteNiftiWrapper(single(pdC), xform, fullfile(outDir,['pdC_.nii.gz']));
    %     disp(['calculate the linear HLF R1b and WF R1b and SIC '])
    %     [HLF,Wf,fh,dd] = hlfCompute(T1,pdC,xform,outDir,MagnetFile,[],0,0,Val/1000);
    %
    %
    %     dtiWriteNiftiWrapper(single(HLF), xform, fullfile(outDir,['HLF_l.nii.gz']));
    %     dtiWriteNiftiWrapper(single(Wf), xform, fullfile(outDir,['Wf_l.nii.gz']));
    %     dtiWriteNiftiWrapper(single(fh), xform, fullfile(outDir,['fh_l.nii.gz']));
    %
    %     [Rb, SIC] = AperntBondT1(T1,brainMask,Wf,MagnetFile);
    %     dtiWriteNiftiWrapper(single(SIC), xform, fullfile(outDir,['SIC_l.nii.gz']));
    %     dtiWriteNiftiWrapper(single(Rb), xform, fullfile(outDir,['R1b_l.nii.gz']));

    %%%%%%%%%%%%%%%%%%%%%%%%%
    
end

return









