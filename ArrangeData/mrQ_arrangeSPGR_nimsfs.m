function mrQ = mrQ_arrangeSPGR_nimsfs(mrQ,inputData_spgr)
% 
%  [mrQ] = mrQ_arrangeSPGR_nimsfs(mrQ,inputData_spgr)
%
% This function will arrange spgr data from a series of scans gathered at the
% CNI. It will take a dataDir (mrQ.RawDir) and organize the data/folders
% within it based upon the series numbers (SPGR_seriesNumbers). 
% It will dump all the data to a raw directory
% (mrQ.outDir).
% 
%     INPUTS:
%                  mrQ:     The mrQ structure
%       inputData_spgr:     The SPGR data.
% 
%    OUTPUTS:
%                  mrQ:     The mrQ structure, updated
%
%
% (C) Stanford University, VISTA Lab [2014]
% Edited by Shai Berman and Jonathan Bain, June-02-2015


%% I. Construct the inputData structures for SEIR and SPGR

% Given mrQ.RawDir, the following function will return all the necessary
% paths to required NIfTI files in the mrQ.inputData_spgr.


    if ~notDefined('inputData_spgr')
    mrQ=mrQ_Set(mrQ,'inputdata_spgr',inputData_spgr);
    end
    %   Go to the nifti hdr, if needed
    if ~isfield(mrQ,'inputdata_spgr') 
        mrQ = mrQ_initInputData_ver2(mrQ,true);
    end
 
 
    


    
%% III. Arrange SPGR data
if ~isfield(mrQ,'MakeNewSPGRRDir');
    mrQ.MakeNewSPGRRDir=0;
end

if (mrQ.MakeNewSPGRRDir==1)
    
    % Initialize counters
    num = 1; ex = 0;
    
    % Make a SPGRDir directory. 
    % If there is one already, this section will make another.
    while ex == 0
        
        SPGRDir=fullfile(mrQ.outDir,['SPGR_' num2str(num)]);
        
        if ( ~exist(SPGRDir,'dir') )
            mkdir(SPGRDir);
            
            ex=1;
            mrQ.MakeNewSPGRRDir=0;
        end
        num=num+1;
    end
    
    mrQ.SPGR = SPGRDir;

    
end

%% IV. Save and return

% This is very important, as it will trick the code into no sorting.
mrQ.ArrangeSPGR_Date = date;

% Save:
save(mrQ.name,'mrQ');

return

