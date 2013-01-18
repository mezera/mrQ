function mrQ = mrQ_Create(path)
%  mrQ = AFQ_Create(path) 
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
% 
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

mrQ.RawDir = path;

% remove spaces and upper case
mrQ.name = fullfile(mrQ.RawDir,'mrQ_params.mat');

% Check for and load mrQ file if it exists
if exist(mrQ.name,'file')
    load (mrQ.name)
    setDefault = 0;
else
    setDefault = 1;
end
    

%% set defult parameters

if setDefault 
    % intiate
    mrQ.clobber = false;

    [~, mrQ.sub] = fileparts(mrQ.RawDir);


    mrQ.arrangeRawFlag=1;
   
    mrQ.channels=[];



% SEIR alignment and fit

    mrQ.complexFlag=0;
    % for comlex data that we not use for now.



    mrQ.useAbs=0;
    % for comlex data that we like to use as magnitude data use 1




% SPGR alignment

    mrQ.coilWeights=1;


    mrQ.alignFlag=1;
    % this to align the SEIR data.we like that unless  if phantom dead or
    % maybe complex (zero)
    mrQ.interp=[];
    mrQ.mmPerVox =[];

    mrQ.refIm=[];
    mrQ.skip=[];
    mrQ.cheack=0;

% fit SPGR

    mrQ.lsq=1; %see for detail inside we recomand the lsq version

    mrQ.runfreesurfer=0;
    mrQ.proclass=0;
    
% A cell array with the SEIR series numbers
mrQ.SEIR_seriesNumbers = {};

% A cell array with the SPGR series numbers
mrQ.SPGR_seriesNumbers = {};


%
save(mrQ.name,'mrQ');

end
