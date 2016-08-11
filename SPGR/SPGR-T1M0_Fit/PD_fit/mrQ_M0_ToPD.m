function [mrQ]=mrQ_M0_ToPD(mrQ)
% function [mrQ]=mrQ_M0_ToPD(mrQ)
%
% This function performs the M0-PD fit. It can be done using T1 single-coil
% data, multi-coil data, or both. The local fits are joined together into
% one PD image.
%
% This function accepts a mrQ structure as its INPUT, and returns an
% updated mrQ structure as its OUTPUT.
%
% See also: mrQ_fitM0boxesCall_T1PD
%           mrQ_fitM0boxesCall
%           mrQ_multicoilM0
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%

%% I. Define inputs
% What method of fit we will use? T1? Multi-coil data? Both?

if ~isfield(mrQ,'PDfit_Method');
    mrQ.PDfit_Method=1;
    % 1. Use only T1 single-coil data to fit (default)
    % 2. Use only multi-coils to fit
    % 3. Use both T1 single-coil data *and* multi-coils to fit
end

if mrQ.PDfit_Method==2
    warning('This fit method was tested primarily with 32 coil channels. Consider using other methods or visit mrQ v1');
end

if mrQ.PDfit_Method~=1
    
    if isfield(mrQ,'calM0_done');
    else
        mrQ.calM0_done=0;
    end
    
    
    if (mrQ.calM0_done==0);
        fprintf('\n Calculate M0 for each coil               \n');
        
        % If we use multi-coil data, we need to calculate the M0 of each coil,
        % build a multi-coil M0 image for the coil's raw data, and then fit T1.
        
        [mrQ.M0combineFile, mrQ.MultiCoilSummedFiles] = mrQ_multicoilM0(mrQ.spgr_initDir,[],[],mrQ.SPGR_niiFile,mrQ.SPGR_niiFile_FA,mrQ);
        mrQ.calM0_done=1;
        save(mrQ.name,'mrQ');
    else
        fprintf('\n Load M0 of each coil               \n');
        
    end
    
end

%% II. Set parameters for the fit

[T1file, M0file,BMfile] = mrQ_get_T1M0_files(mrQ);
if mrQ.PDfit_Method~=1
    M0file=mrQ.M0combineFile;
end

[mrQ.opt_logname] = mrQ_PD_Fit_saveParams(mrQ.spgr_initDir,mrQ.sub,mrQ.PolyDeg,M0file,T1file,BMfile,mrQ.PDfit_Method,mrQ);
save(mrQ.name,'mrQ');

%% III. Perform the fit

if ~isfield(mrQ,'SPGR_PDfit_done');
    mrQ.SPGR_PDfit_done=0;
end

if mrQ.SPGR_PDfit_done==0
if mrQ.PDfit_Method==1
    % Note that we fit each local M0 area (with T1 input). This is fast.
    % We don't need the grid, but we can expedite it with parloop.
    fprintf('\n fitting M0 ...              \n');
    
    mrQ_fitM0boxesCall_T1PD(mrQ.opt_logname);
    mrQ.SPGR_PDfit_done=1;
    save(mrQ.name,'mrQ');
    
elseif mrQ.PDfit_Method==2 ||  mrQ.PDfit_Method==3
    % Multi-coil fit.
    % This will be better with parallel computing like SunGrid.
    % ** We need to test parloop here also. **
    fprintf('\n fitting multi coil M0 ...              \n');
    
    mrQ_fitM0boxesCall(mrQ.opt_logname,mrQ.SunGrid);
    mrQ.SPGR_PDfit_done=1;
    save(mrQ.name,'mrQ');
    
end
else
            fprintf('\n using previously fitted PD boxes               \n');
end

%% IV. Build the local fits
%Join the local overlap area to one PD image.
fprintf('\n Building the M0 map...              \n');

if ~isfield(mrQ,'RepErrThreshold')
   mrQ.RepErrThreshold=[];
end
if ~isfield(mrQ,'PrcCutOff')
   mrQ.PrcCutOff=[];
end
if ~isfield(mrQ,'ErrorThresh')
   mrQ.ErrorThresh=[];
end

                              
mrQ.opt_logname=mrQ_buildPD(mrQ.opt_logname,[],[],mrQ.RepErrThreshold,mrQ.PrcCutOff,mrQ.ErrorThresh);
%

