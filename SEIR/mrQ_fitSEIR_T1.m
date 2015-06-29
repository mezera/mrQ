function [mrQ]= mrQ_fitSEIR_T1(SEIRdir,outDir,checkData,mrQ)
% 
%  [mrQ]=mrQ_fitSEIR_T1(SEIRdir,outDir,[checkData=1])
% 
% Loads SEIR_Dat.mat file from the directory 'data' in the SEIRdir, and
% performs T1 mapping, also displays the results for the user to check if
% checkData = 1.
% 
% 
% INPUTS:
%   SEIRdir -     Path to your desired T1 image. Use the same path as in
%                 getSEIR.
% 
%   checkData -   If you want to visually check the data leave empty or set
%                 to 1. To not check data set to 0.
%   mrQ      -     information structure
% 
% OUTPUTS:
%     mrQ      -    information structure, updated
%   
%   This function will save the fitted images in a new directory called
%   "fitT1_GS" under the SEIRdir. the name of the saved files in  the
%   output strctures
%
% WEB RESOURCES:
%   web('http://white.stanford.edu/newlm/index.php/Quantitative_Imaging','-browser');
% 
% 
% EXAMPLE USAGE:
%   SEIRdir = '/biac2/wandell2/data/WMDevo/adult/109_AL/QuantitativeImaging/20110622_0582/SEIR_epi_1'
%   mrQ_fitSEIR_T1(SEIRdir,[],[],[],mrQ);
%   
 
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University
% 

%% Check INPUTS

% This makes many assumptions about the paths. This should be flexible. All
% we need is the data - the input should just point to that. The function
% can then move one directory up from there and create the fit* directory.

if notDefined('SEIRdir') || ~exist(SEIRdir,'dir')
    SEIRdir = uigetdir(pwd,'Select your SEIR base data directory');
end

seirDataDir = fullfile(SEIRdir,'data');

if notDefined('outDir')
    outDir = fullfile(SEIRdir,'fitT1_GS');
end

if ~exist(outDir,'dir'), mkdir(outDir); end


if notDefined('checkData') || isempty(checkData) || checkData > 1
    checkData = 1;
end

% Which algorithm to use: real 
% nonlinear least squares polarity restoration. 

method = 'NLSPR'; % magnitude data

close all


%% Load the data file (this will load 'data', 'xform', and 'extra')

% We should now be pointed to this file, and not to the directory. It is
% much easier to work with filenames than with directories.
filename = 'SEIR_Dat';
loadStr = fullfile(seirDataDir, filename);

outName = ['T1Fit' method '_' filename];  
saveStr = fullfile(outDir, outName);


%% Perform T1 fits 

% Estimates T1 together with:
%   NLSPR: a and b parameters to fit the data to |a + b*exp(-TI/T1)|
T1FitExperimentData(loadStr, saveStr, method, checkData);


% Load the fitted data 
% ** We have to do it this way because T1FitExperimentData does not return
% the data as desired.
load(saveStr)
load(loadStr)


%% Save Nifti

% The 'll_T1' variable comes from the function T1FitExperimentData - it
% could be returned from there instead of loading the file containing it. 

T1file=[saveStr '_T1.nii.gz'];
resnormfile=[saveStr '_T1FitResid.nii.gz'];
dtiWriteNiftiWrapper(single(ll_T1(:,:,:,1)), xform, T1file); %#ok<NODEF>
dtiWriteNiftiWrapper(single(ll_T1(:,:,:,4)), xform, resnormfile);

mrQ.SEIR_epi_T1file=T1file;
mrQ.SEIR_epi_resnormfile=resnormfile;
mrQ.SEIR_epi_fitFile=saveStr;


return
