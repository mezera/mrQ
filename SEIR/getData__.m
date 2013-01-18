% getData.mat
%
% Generates the .mat from the DICOM images. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    

clear all
close all

%T1path = '/biac3/wandell5/data/relaxometry/KK02232010/';%'../';
%T1path = '/biac3/wandell5/data/relaxometry/2010_06_29_rfd_T1_mapping/';
%T1path = '/biac3/wandell5/data/relaxometry/amGS/';
%T1path = '/biac3/wandell5/data/relaxometry/AKGS/';
%T1path = '/biac3/wandell5/data/relaxometry/am_GS_3T_21397_11022010/';
%T1path = '/biac3/wandell5/data/relaxometry/AM3TGS12082010/';
% T1path = '/biac3/wandell5/data/relaxometry/RB3TGS12082010/';
% T1path = '/biac3/wandell5/data/relaxometry/phantom11062010/5461_11052010/';
T1path = '/biac3/wandell5/data/relaxometry/RM20110531/';

loadpath = [T1path 'data/'];%'singleslicedicomfiles/']; 
savename = [T1path 'data/' 'SEIR3T_Fit']; 

currpath = pwd;

if(~exist(loadpath,'dir')), mkdir(loadpath); end



cd(loadpath)
addpath([T1path '/mfiles']);

d = dicomLoadAllSeries('.');
cd(currpath)

complexFlag = 0;%1; % 1: complex data; 0: magnitude data

nbrow = size(d(1).imData,1);
nbcol = size(d(1).imData,2);
nbslice = size(d(1).imData, 3);

nbseries = length(d);

data = zeros(nbrow,nbcol,nbslice,nbseries); % Complex data
extra.tVec = zeros(1,nbseries); % One series corresponds to one TI

for k = 1:nbseries
	dataTmp = d(k).imData;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice
		if (complexFlag)
			data(:,:,ss,k) = dataTmp(:,:,3+(ss-1)*4)+i*dataTmp(:,:,4+(ss-1)*4);% Complex
		else
			%data(:,:,ss,k) = dataTmp(:,:,1+(ss-1)*4); % Magnitude A.M
			%7/8/2010: I think it is wrong and not getting the right slices
			%so i chance it. i don't get the reson of the original code.
            data(:,:,ss,k) = dataTmp(:,:,ss); % Magnitude
        end      
	end
	extra.tVec(k) = d(k).inversionTime;
end 

% This is where you would correct the phase of certain datapoints (if
% Prescan had been used between scans)
% e.g., data(:,:,3) = -data(:,:,3); 
%    or data(:,:,:,3) = -data(:,:,:,3); if you have more than one slice
% !!! (:,:,1) corresponds to the first series acquired, not to the smallest TI !!!

extra.T1Vec = 1:5000; % This can be reduced to speed up the code

TI = extra.tVec;
save(savename,'data','extra')
