% mainCurr.m
%
% Loads .mat file from the directory 'data' in the T1path, 
% performs T1 mapping, and displays the results 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

%clear all
close all

%T1path = '/biac3/wandell5/data/relaxometry/KK02232010/';%'../';
%T1path = '/biac3/wandell5/data/relaxometry/2010_06_29_rfd_T1_mapping/';
%T1path = '/biac3/wandell5/data/relaxometry/amGS/';
%T1path = '/biac3/wandell5/data/relaxometry/AKGS/';
%T1path = '/biac3/wandell5/data/relaxometry/am_GS_3T_21397_11022010/';
%T1path = '/biac3/wandell5/data/relaxometry/AM3TGS12082010/';
% T1path = '/biac3/wandell5/data/relaxometry/RB3TGS12082010/';
 %T1path = '/biac3/wandell5/data/relaxometry/phantom11062010/5461_11052010/';
T1path = '/biac3/wandell5/data/relaxometry/RM20110531/';
% Where to find the data
loadpath = [T1path 'data/'];

%filename = 'SE1p5T_singleslice';
%filename = 'SE1p5T_ThreeSlice_Mg';
filename = 'SEIR3T_Fit';

complexFlag =0; % 1: complex data; 0: magnitude data

% Which algorithm to use
if (complexFlag)
	method = 'NLS'; % complex data
else
	method = 'NLSPR'; % magnitude data
end

% Where to save the data
%savepath = [T1path 'fitdata/'] 
savepath = [T1path 'fitT1_GS/'];

if(~exist(savepath,'dir')), mkdir(savepath); end
loadStr = [loadpath filename]
saveStr = [savepath 'T1Fit' method '_' filename]

% Perform fit
T1FitExperimentData(loadStr, saveStr, method);

% Display results
T1FitDisplayBrain(loadStr, saveStr, method);
%T1FitDisplayPhantom(loadStr, saveStr, method);

