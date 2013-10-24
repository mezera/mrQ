function T1FitExperimentData(loadStr, saveStr, method, checkData)
% 
% T1FitExperimentData(loadStr, saveStr, method,[checkData=1])
% 
% loadStr:    Data file to be loaded
% saveStr:    Where results of the fit will be saved along with filename
% method:     What fitting method to use: NLS or NLSPR
% checkData:  If you want to visually check the data leave empty or set to
%             1. To not check data set to 0. 
%
% Estimates T1 together with:
%   NLS:   a and b parameters to fit the data to a + b*exp(-TI/T1)
%   NLSPR: a and b parameters to fit the data to |a + b*exp(-TI/T1)
%
%
%  (c) Board of Trustees, Leland Stanford Junior University

% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009


%% Check INPUTS

if notDefined('loadStr') || ~exist([loadStr '.mat'],'file')
    loadStr = uigetdir(pwd,'Select the directory containing your SEIR data');
end

if notDefined('checkData') || isempty(checkData) || checkData > 1
    checkData = 1;
end

% Hardcode the saving of the fit data. We could return what we need and
% make this an option. 
savefitdata = 1; 

 
%% Load the data and setup the commands that will perform the fit 

% (see getSEIR.m mrQ_getSEIR.m) This file (SEIR_Dat.mat) file contains
% all the DICOM image data in the variable 'data' (4-D double) as well as
% the transform used to align the SEIR series in the variable 'xform'. Also
% in the .mat file is the variable 'extra' which contains 'tVec' - (1xN)
% the inversion times for each series (N) and 'T1Vec' - [1:5000].
load(loadStr);

switch method
  case {'NLS','NLSPR'}
    nlsS = getNLSStruct(extra,1); %#ok<NASGU>
end

switch method
    case {'NLS'}
        nbrOfFitParams = 4; % Number of output arguments for the fit
        fitStr   = ['[T1Est bEst aEst res] = T1Fit' method '(data(jj,:),nlsS);'];
        storeStr = 'll_T1(jj,:) = [T1Est bEst aEst res];';
        clearStr = 'clear jj T1Est bEst aEst res';
    case {'NLSPR'}
        % We're only working with the magnitude data
        data = abs(data); %#ok<NODEF>
        nbrOfFitParams = 4; % Number of output arguments for the fit
        fitStr   = ['[T1Est bMagEst aMagEst res] = T1Fit' method '(data(jj,:),nlsS);'];
        storeStr = 'll_T1(jj,:) = [T1Est bMagEst aMagEst res];';
        clearStr = 'clear jj T1Est bMagEst aMagEst res';
end

dims  = size(data);
nbrow = size(data,1);
nbcol = size(data,2);

if numel(dims) > 3
    nbslice = dims(3); % Check number of slices 
else
    nbslice = 1;
    tmpData(:,:,1,:) = data; % Make data a 4-D array regardless of number of slices
    data = tmpData;
    clear tmpData;
end


%% Check Data

% Adding the option to check the data manually. Manual data checking is not
% always practical to do given that many datasets may be processed in
% sequence.

% This cell may have some errors. It has been noted that we're getting
% errors if I choose slice 12 - or if I say that the mask is no good. *
% Need to explain what we're doing here and how to choose a slice, and how
% to evaluate the graphs that appear *

if checkData
    % The data needs to lie along an exponential.
    % Refer to getData.m if that is not the case.
    while (true)
        if (nbslice ~= 1)
            prompt = sprintf('For which slice would you like to check the data?\nEnter number 1 to %s [0 for no check]. \n',num2str(nbslice));
            zz = inputdlg(prompt, 'Data Check');
            zz = cast(str2double(zz), 'int16');
            if (isinteger(zz) && zz >= 0 && zz <= nbslice), break;
            end
        else
            zz = inputdlg('Enter 1 to check the data, 0 for no check.', 'Data check');
            zz = cast(str2double(zz), 'int16');
            if (isinteger(zz) && zz >= 0 && zz <= nbslice), break;
            end
        end
    end
    
    if (zz ~=0)
        sliceData = squeeze(data(:,:,zz,:));
        msgbox('Click on one point to check that the time series looks like an exponential. CTRL-click or right-click when done. Refer to getData.m if the data do not look correct.' , 'Data Check Instuctions','help','modal');
        TI = extra.tVec;
        plotData(real(sliceData),TI);
    end
    
    close all
end

%%

dataOriginal = data;


%% Mask the background

% For each slice we mask points dimmer than maskFactor*the brightest point. 
mesures=size(data,4);
maskFactor = 0.03;
mask = zeros(nbrow, nbcol, nbslice,mesures);

%[u,v] = max(extra.tVec); %#ok<ASGLU>
% The intensity mask is taken with respect to the image with the largest TI, 
% where we expect the magnetization to have almost fully recovered. 
for v=1:size(data,4)
for kk = 1:nbslice
	maskTmp = mask(:,:,kk);
	%maskTmp = medfilt2(maskTmp); % remove salt and pepper noise
	maskThreshold = maskFactor * max(max(abs(data(:,:,kk,v))));
    maskTmp=logical(abs(data(:,:,kk,v))>maskThreshold);
	maskTmp(find(abs(data(:,:,kk,v)) > maskThreshold)) = 1; %#ok<FNDSB>
	mask(:,:,kk,v) = maskTmp;
	clear maskTmp
end
end
mask=(sum(mask,4));
mask=logical(mask==mesures);
%let fit also around the mask so we won't have holes (that can help registration to the SPGR)
D3=size(mask,3);
if D3>1
mask1=logical(smooth3(mask));%,'box',[5 5 5]
else
    mask1=mask;
end
   

    
 for i=1:size(mask1,3)
 mask1(:,:,i)=imfill(mask1(:,:,i),'holes');
 end;
 mask=double(mask1);
% but exclude zeros nan and inf voxels
 for i=1:length(extra.tVec);
     tmp=data(:,:,:,i);
     mask(isnan(tmp))=0;
     mask(isinf(tmp))=0;
      mask(tmp<0)=0;
 end


maskInds = find(mask);    


%% Check Data: Mask inspection

if checkData
    % Show the mask to the user and ask them to ok it.
    showMontage(mask) ;title('Mask')
    % figure,
    % imshow(abs(mask(:,:,1)),[])
    % title('Mask for slice 1')
    customFormat; % Set custom graph properties
    
    if (input('Enter 1 if mask is appropriate, 0 if not --- ')==0)
        error('Change maskFactor in T1FitExperimentData.m')
    end
end


%% Do the actual fitting (see fitStr above)

% How many voxels to process before printing out status data
nVoxAll            = length(maskInds);
numVoxelsPerUpdate = min(floor(nVoxAll/10),1000); 
						   
ll_T1  = zeros(nVoxAll, nbrOfFitParams);
nSteps = ceil(nVoxAll/numVoxelsPerUpdate); % Number of status reports

for ii=1:size(data,4)
    tmpVol = data(:,:,:,ii);
    tmpData(:, ii) = tmpVol(maskInds)';
end

clear tmpVol;
data = tmpData; %#ok<NASGU>
clear tmpData;

startTime = cputime;
fprintf('Processing %d voxels.\n', nVoxAll);
h = waitbar(0, sprintf('Processing %d voxels', nVoxAll)); 

for ii=1:nSteps
  curInd = (ii-1)*numVoxelsPerUpdate+1;
  endInd = min(curInd+numVoxelsPerUpdate,nVoxAll);
  
  for jj = curInd:endInd
    
    eval(fitStr)    % Do the fit
    
    eval(storeStr); % Store the data
    
  end
  
  waitbar(ii/nSteps, h, sprintf('Processing %d voxels, %g percent done...\n',nVoxAll,round(endInd/nVoxAll*100)));
end

% Clean up
eval(clearStr)
close(h);
timeTaken = round(cputime - startTime);
fprintf('Processed %d voxels in %g seconds.\n',nVoxAll, timeTaken);

dims = [size(mask) 4]; %#ok<NASGU>
im   = zeros(size(mask));

for ii=1:nbrOfFitParams
    im(maskInds) = ll_T1(:,ii);
    T1(:,:,:,ii) = im; %#ok<AGROW>
end

% Going back from a numVoxels x 4 array to nbrow x nbcol x nbslice
ll_T1 = T1;


%% Store ll_T1 and mask in saveStr
% ll_T1 has four parameters for each voxel, namely:
% (1) T1 
% (2) 'b' parameter from the model 
% (3) 'a' parameter
% (4) residual from the fit

if (savefitdata)
	save(saveStr,'ll_T1','mask','nlsS')
end


%% Data check: Fit inspection

if checkData
    % Quick check of the fit
    TI    = extra.tVec;
    nbtp  = 20;
    timef = linspace(min(TI),max(TI),nbtp);
    
    % Inserting a short pause, otherwise some computers seem to have problems
    pause(1);
    
    zz = 0;
    while(true)
        if (nbslice ~= 1)
            disp('For which slice would you like to check the fit?');
            zz = input(['Enter number 1 to ' num2str(nbslice) '.  0 for no check --- '], 's');
            zz = cast(str2double(zz), 'int16');
            if (isinteger(zz) && zz >= 0 && zz <= nbslice)
                break;
            end
        else
            zz = input('Enter 1 to check the fit, 0 for no check --- ', 's');
            zz = cast(str2double(zz), 'int16');
            if (isinteger(zz) && zz >= 0 && zz <= nbslice)
                break;
            end
        end
    end
    
    if (zz ~= 0)
        sliceData = squeeze(dataOriginal(:,:,zz,:));
        datafit = zeros(nbrow,nbcol,nbtp);
        switch method
            case{'NLS'}
                for kk = 1:20
                    datafit(:,:,kk) = ll_T1(:,:,zz,3) + ...
                        ll_T1(:,:,zz,2).*exp(-timef(kk)./ll_T1(:,:,zz,1));
                end
            case{'NLSPR'}
                for kk = 1:20
                    datafit(:,:,kk) = abs(ll_T1(:,:,zz,3) + ...
                        ll_T1(:,:,zz,2).*exp(-timef(kk)./ll_T1(:,:,zz,1)));
                end
        end
        
        % Check Data: Fit inspection
        disp('Click on one point to check the fit. CTRL-click or right-click when done')
        plotData(real(sliceData),TI,real(datafit),ll_T1);
        close all
    end
end
