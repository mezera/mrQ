function mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir)
% 
%  [mrQ] = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir)
%
% Arrange data from a series of scans gathered at the CNI. This function
% will take a dataDir (mrQ.RawDir) and organize the data/folders within it
% based upon the series numbers (SEIR_seriesNumbers & SPGR_seriesNumbers).
% It will dump all the data to a raw directory (mrQ.outDir).
% 
%     INPUTS:
%      
%       inputData_spgr:     The SPGR data.
% 
%       inputData_seir:     The SEIR data.
% 
% 
% (C) Stanford Univsersity, VISTA Lab [2014]
% Edited by Shai Berman and Jonathan Bain, June-02-2015


%% Construct the inputData structures for SEIR and SPGR

% Given mrQ.RawDir, the following function will return all the necessary
% paths to required nifti files in the mrQ.inputData_spgr and
% mrQ.inputData_seir structures.
if ~notDefined('inputData_seir')
    mrQ=mrQ_Set(mrQ,'inputdata_seir',inputData_seir);
end
if ~notDefined('inputData_spgr')
    mrQ=mrQ_Set(mrQ,'inputdata_spgr',inputData_spgr);
end
%   go to the nifti hdr if needed
if ~isfield(mrQ,'inputdata_spgr') || ~isfield(mrQ,'inputdata_seir')
    mrQ = mrQ_initInputData(mrQ);
end

%% Arrange SEIR data

if ~isfield(mrQ,'MakeNewSEIRDir');
    mrQ.MakeNewSEIRDir=0;
end

if (mrQ.MakeNewSEIRDir==1)
    
    % Initialize counters
    num = 1; ex = 0;
    
    % Make a SEIR_epi directory. 
    % If there is one already, this section will make another.
    while ex == 0
        SEIRepiDir = fullfile(mrQ.outDir,['SEIR_epi_' num2str(num)]);
        
        if ( ~exist(SEIRepiDir,'dir') )
            mkdir(SEIRepiDir);
            
            SEIRepiDir_dat = fullfile(SEIRepiDir,'data');
            mkdir(SEIRepiDir_dat);
            
            SEIRepiDir_fit = fullfile(SEIRepiDir,'fitT1_GS');
            mkdir(SEIRepiDir_fit);
            ex=1;
            mrQ.MakeNewSEIRDir=0;
            
        end
        num=num+1;
    end
    
    mrQ.SEIRepiDir=SEIRepiDir;
%     mrQ.SEIR_done=0;

end 
    
%% Arrange SPGR data
if ~isfield(mrQ,'MakeNewSPGRRDir');
    mrQ.MakeNewSPGRRDir=0;
end

if (mrQ.MakeNewSPGRRDir==1)
    
    % Initialize counters
    num = 1; ex = 0;
    
    % Make a SPGRDir Dir. If there is one already it will make another.
    while ex == 0
        
        SPGRDir=fullfile(mrQ.outDir,['SPGR_' num2str(num)]);
        
        if ( ~exist(SPGRDir,'dir') )
            mkdir(SPGRDir);
            
            SPGRDir_dat = fullfile(SPGRDir,'data'); mkdir(SPGRDir_dat);
            ex=1;
            mrQ.MakeNewSPGRRDir=0;
        end
        num=num+1;
    end
    
    mrQ.SPGR = SPGRDir;
%     mrQ.SPGR_T1fit = 0;
%     mrQ.SPGR_init = 0;
%     mrQ.SPGR_PDfit = 0;
%     mrQ.SPGR_PDBuild = 0;
    
end

%% Save and return

% This is very important, as it will trick the code into no sorting.
mrQ.Arange_Date = date;

% Save:
save(mrQ.name,'mrQ');

return

