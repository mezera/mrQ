function [mrQ]=mrQ_M0_ToPD(mrQ)

%1. Define inputs



%% What method of fit we will use? Will we use multi coil data?


if ~isfield(mrQ,'PDfit_Method');
    
    mrQ.PDfit_Method=1; % 1. will use only T1 to fit (default old 5); 2. Will use only multi coils to fit (old 2); 3. Will use both multi coils and T1 to fit (old 1);
end

if mrQ.PDfit_Method~=1
    
    if isfield(mrQ,'calM0_done');
    else
        mrQ.calM0_done=0;
    end
    
    if (mrQ.calM0_done==0);
        fprintf('\n calculate M0 for each coil               \n');
        
        %  if we using multi coil data we need to calculate the M0 of each coil
        %build a multi coil M0 image for the coils raw data and then fitted T1
        [mrQ.M0combineFile] = mrQ_multicoilM0(mrQ.spgr_initDir,[],[],mrQ.SPGR_niiFile,mrQ.SPGR_niiFile_FA,mrQ);
        
        save(mrQ.name,'mrQ');
    else
        fprintf('\n load M0 of each coil               \n');
        
    end
end

%%
%% set parameters for the fit

[T1file, M0file,BMfile] = mrQ_get_T1M0_files(mrQ);

[mrQ.opt_logname] = mrQ_PD_Fit_saveParams(mrQ.spgr_initDir,mrQ.sub,mrQ.PolyDeg,M0file,T1file,BMfile,mrQ.PDfit_Method,mrQ);
save(mrQ.name,'mrQ');

%

%%
%3 fit

if ~isfield(mrQ,'SPGR_PDfit_done');
    mrQ.SPGR_PDfit_done=0;
end


%%

if mrQ.PDfit_Method==1
  %Note fit each local M0 area (with T1 input). this is fast. We don't need the
% grid we can expedite it with parloop  
    mrQ_fitM0boxesCall_T1PD(mrQ.opt_logname);
elseif mrQ.PDfit_Method==2 ||  mrQ.PDfit_Method==3
    % multi coil fit, this will be better with parallel computing like
    % sunGrid. We need to test Parloop here also.
            mrQ_fitM0boxesCall(mrQ.opt_logname,mrQ.SunGrid);
    
end
%%


    

%% Build the local fits
%Join the local overlap area to one PD image.

mrQ.opt=mrQ_buildPD_ver2(mrQ.opt_logname,[],[],[],[],0.01);
%

%4 build
