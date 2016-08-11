function mrQ = mrQ_arrangeSEIR_nimsfs(mrQ,inputData_seir)
% 
%  [mrQ] = mrQ_arrangeseir_nimsfs(mrQ,inputData_spgr,inputData_seir)
%
% This function will arrange data from a series of scans gathered at the
% CNI. It will take a dataDir (mrQ.RawDir) and organize the data/folders
% within it based upon the series numbers (SEIR_seriesNumbers). 
% It will dump all the data to a raw directory
% (mrQ.outDir).
% 
%     INPUTS:
%                  mrQ:     The mrQ structure
%       inputData_seir:     The SEIR data.
% 
%    OUTPUTS:
%                  mrQ:     The mrQ structure, updated
%
%
% (C) Stanford University, VISTA Lab [2014]
% Edited by Shai Berman and Jonathan Bain, June-02-2015


%% I. Construct the inputData structures for SEIR and SPGR

% Given mrQ.RawDir, the following function will return all the necessary
% paths to required NIfTI files in the
% mrQ.inputData_seir structures.



   if ~notDefined('inputData_seir')
    mrQ=mrQ_Set(mrQ,'inputdata_seir',inputData_seir);
    end
    %   Go to the nifti hdr, if needed
    if ~isfield(mrQ,'inputdata_seir') 
        mrQ = mrQ_initInputData_ver2(mrQ,false,true);
    end





%% IV. Save and return

% This is very important, as it will trick the code into no sorting.
mrQ.ArrangeSEIR_Date = date;

% Save:
save(mrQ.name,'mrQ');

return

