function mrQ = mrQ_Create(path,name,outDir)
%  mrQ = mrQ_Create(path)
%
% This function will load data into a mrQ structure and set default
% parameters. if the mrQ file exsists withing the directory where this
% funciton is called, then it will load it
%
%
% INPUTS:
%     path -  The path is a string of location wehre the structure is saved.
%             the path will be also set as the rawDir. Therefore it is
%             usful to make path the location of the raw image.
%     name -  a string that will be saved as the name of the mrQ structure. The defult
%             is "mrQ_params"

%
% EXAMPLE USAGE:
%
% mrQ = mrQ_Create('/biac2/wandell2/data/WMDevo/ambliopia/sub7/QuantativeImaging/20121102_3488')
%
% (c) STANFORD UNIVERSITY, VISTA LAB
%

%%

if notDefined('path') || ~exist(path,'dir')
    path = uigetdir(pwd,'Choose mrQ directory');
end

if notDefined('outDir') || ~exist(outDir,'dir')
    outDir = path;
else
    mrQ.outDir = outDir;
end

mrQ.RawDir = path;
mrQ.outDir = outDir;

if notDefined('name')
    mrQ.name = fullfile(mrQ.outDir,'mrQ_params.mat');
else
    mrQ.name = fullfile(mrQ.outDir,[name '.mat']); 
end

% Check for and load mrQ file if it exists
if exist(mrQ.name,'file')
    load (mrQ.name)
    setDefault = 0;
    fprintf('mrQ file found:\n\tLoading %s.\n',mrQ.name);
else
    setDefault = 1;
end


%% set defult parameters

if setDefault
    %% intiate
    mrQ.clobber = false;
    [~, mrQ.sub] = fileparts(mrQ.RawDir);
    mrQ.arrangeRawFlag=1;
    mrQ.channels=[];
    %% SEIR
    % SEIR alignment and fit
    mrQ.MakeNewSEIRDir=1;
    
    mrQ.SEIR_done=0;
    mrQ.complexFlag=0;
    % for comlex data that we not use for now.
    mrQ.alignFlag=1;
    mrQ.useAbs=0;
    % for comlex data that we like to use as magnitude data use 1
    % A cell array with the SEIR series numbers
    mrQ.SEIR_seriesNumbers = {};
    
    %% SPGR
    mrQ.MakeNewSPGRRDir=1;
    
    mrQ.SPGR_init_done=0;
    mrQ.SPGR_coilWeight_done=0;
    mrQ.SPGR_T1fit_done=0;
    
    % SPGR alignmen
    mrQ.coilWeights=1;
    
    % this to align the SEIR data.we like that unless  if phantom dead or
    % maybe complex (zero)
    mrQ.interp=[];
    mrQ.mmPerVox =[];
    
    mrQ.refIm=[];
    mrQ.skip=[];
    mrQ.permution=0;
    % A cell array with the SPGR series numbers
    mrQ.SPGR_seriesNumbers = {};
    
    % brack and cheaks
    mrQ.cheack=0;
    mrQ.brakeAfterVisualization=0;
    mrQ.brakeAfterT1=0;
    mrQ.viewbrake=0;
    %%  T1 fit of SPGR
    mrQ.lsq=1; %see for detail inside we recomand the lsq version
    
    %% Segmentation
    mrQ.segmentaion=0;
    mrQ.runfreesurfer=0;
    
    
    
    %% fit PD
    mrQ.calM0_done=0;
    SPGR_PDfit_done=0;
    mrQ.SPGR_PDBuild_done=0;
    
    mrQ.PolyDeg=3;
    mrQ.SunGrid=1;
    mrQ.proclus=1;
    %
    save(mrQ.name,'mrQ');
    
end
