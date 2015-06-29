function mrQ = mrQ_Create(path,name,outDir)
%  mrQ = mrQ_Create(path)
%
% This function will load data into a mrQ structure and set default
% parameters. If the mrQ file exists within the directory where this
% function is called, then the function will load the file.
%
%
% INPUTS:
%     path -  The path is a string of location where the structure is saved.
%             The path will be also set as the rawDir. Therefore it is
%             useful to make the path the location of the raw image.
%     name -  A string that will be saved as the name of the mrQ structure. 
%             The default name is "mrQ_params".
%     outDir- 
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
    [~, mrQ.sub] =fileparts(tempname);
    
  
    %% SEIR
    % SEIR alignment and fit
    mrQ.MakeNewSEIRDir=true;
        mrQ.alignFlag=true;

    mrQ.SEIR_done=false;
  %  mrQ.complexFlag=false;
    % for complex data that we do not use for now.
   % mrQ.useAbs=false;
    % for complex data that we'd like to use as magnitude data. Use 1
   
    % Create a cell array to hold the SEIR series numbers:
    mrQ.SEIR_seriesNumbers = {};
    
    %% SPGR
    mrQ.MakeNewSPGRRDir=true;
    
    mrQ.SPGR_init_done=false;
    
%    mrQ.SPGR_coilWeight_done=0;

    mrQ.SPGR_T1fit_done=false;
    
    % SPGR alignment
   % mrQ.coilWeights=1;
    
    % This is to align the SEIR data.
    %we like that unless  if phantom dead or maybe complex (zero)
    mrQ.interp=[];
    mrQ.mmPerVox =[];
    
    mrQ.refIm=[];
    mrQ.skip=[];
    
    mrQ.permution=false;
    % Create a cell array for the SPGR series numbers:
    mrQ.SPGR_seriesNumbers = {};
    
    %% for debugging
    % brake and checks
%     mrQ.check=0;
%     mrQ.brakeAfterVisualization=0;
%     mrQ.brakeAfterT1=0;
%     mrQ.viewbrake=0;

    %%  T1 fit of SPGR
    mrQ.lsq=false; %see for detail inside we recomand the lsq version
    mrQ.LW=true;  %or the WL fits
    
    %% Segmentation
    mrQ.segmentation=false;
    mrQ.runfreesurfer=false; % the defult is FSL
    
    %% fit PD
    mrQ.calM0_done=false;
    mrQ.SPGR_PDfit_done=false; %   
    mrQ.SPGR_PDBuild_done=false;
    
    mrQ.PolyDeg=3; 
    mrQ.SunGrid=false;
    
   % mrQ.proclus=1;
    %
    save(mrQ.name,'mrQ');
    
end
