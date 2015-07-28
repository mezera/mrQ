function [T1file, M0file,BMfile]=mrQ_get_T1M0_files(mrQ,whT1,whM0,whBM)
% function [T1file, M0file,BMfile]=mrQ_get_T1M0_files(mrQ,whT1,whM0,whBM)
%
% Search for the T1, M0 and Brain Mask files. Can be one of three file
% names, depending on the fit process.
%
%  ~INPUTS~
%          mrQ: The mrQ structure
%         whT1: Checker if the T1 file exists
%         whM0: Checker if the M0 file exists
%         whBM: Checker if the Brain Mask file exists
%
%  ~OUTPUTS~
%       T1file: Location of the T1 file.
%       M0file: Location of the M0 file.
%       BMfile: Location of the Brain Mask file.
%


%% I. Definitions
if notDefined('whT1')
    whT1=true;
end

if notDefined('whM0')
    whM0=true;
end
if notDefined('whBM')
    whBM=true;
end

if whT1
    T1file=[];
    
    
%% II. Looking for the T1 file  
% Look for linear fit, weighted linear, non-linear
% Then, check manually
% If all this fails, present error message.

%% II-a. Is there a linear fit T1 file?
    %Find T1
    
    if isfield(mrQ,'T1_LFit')
        if  exist(mrQ.T1_LFit,'file')
            T1file= mrQ.T1_LFit;
        end
    end
    
    %  the same but with larger mask
    if isfield(mrQ,'T1_LFit_HM')
        if  exist(mrQ.T1_LFit_HM,'file')
            T1file= mrQ.T1_LFit_HM;
        end
    end
    
    %  Maybe after B1 correction?
    if isfield(mrQ,'T1_B1_LFit')
        if exist(mrQ.T1_B1_LFit,'file') 
            T1file= mrQ.T1_B1_LFit;
        end
    end
    
%% II-b. Is there a weighted linear T1 file?
    if isfield(mrQ,'T1_B1_LWFit')
        
        if exist(mrQ.T1_B1_LWFit,'file')
            T1file= mrQ.T1_B1_LWFit;
        end
    end
    
%% II-c. Is there an non-linear T1 file?
    if isfield(mrQ,'T1_B1_lsqFit')
        
        if  exist(mrQ.T1_B1_lsqFit,'file')
            T1file= mrQ.T1_B1_lsqFit;
        end
    end
    
%% II-d. Select manually
    if isempty(T1file)
        T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
    end
    
%% II-e. Still can't find it (error)
    
    if isempty(T1file)
        error('Cannot find the T1file')
    end
else
    T1file=[];
end

%
%%%
%%%%%
%%%
%

%% III. Looking for the M0 file  
% Look for linear fit, weighted linear, non-linear
% Then, check in the outputDir
% If all this fails, present error message.

%% III-a. Is there a linear fit M0 file?
%Find M0
if whM0
    M0file=[];
    if isfield(mrQ,'M0_LFit')
        if  exist(mrQ.M0_LFit,'file') 
            M0file= mrQ.M0_LFit;
        end
    end
    
    %  the same but with larger mask
    if isfield(mrQ,'M0_LFit_HM')
        if  exist(mrQ.M0_LFit_HM,'file') 
            M0file= mrQ.M0_LFit_HM;
        end
    end
    
    
    %  Maybe after B1 correction?
    if isfield(mrQ,'M0_B1_LFit')
        if  exist(mrQ.M0_B1_LFit,'file') 
            M0file= mrQ.M0_B1_LFit;
        end
    end
      
%% III-b. Is there a weighted-linear M0 file?
    
    if isfield(mrQ,'M0_B1_LWFit')
        
        if exist(mrQ.M0_B1_LWFit,'file')
            M0file= mrQ.M0_B1_LWFit;
        end
    end
%% III-c. Is there a non-linear M0 file?
    
    if isfield(mrQ,'M0_B1_lsqFit')
        
        if  exist(mrQ.M0_B1_lsqFit,'file')
            M0file= mrQ.M0_B1_lsqFit;
        end
    end
    
    % Are we using multi coils data. if yes this will be the M0
     if isfield(mrQ,'M0combineFile')
        
        if  exist(mrQ.M0combineFile,'file')
            M0file= mrQ.M0combineFile;
        end
    end

%% III-d. Maybe it is already in the outputDir?
    if isfield(mrQ, 'maps')
        if isfield(mrQ.maps,'T1path')
            if exist(mrQ.maps.T1path,'file')
                T1file= mrQ.maps.T1path;
            end
        end
    end
    
    % Note that we are not moving the M0 image to the output Dir, so we are
    % not looking for it there.
    
    if isempty(M0file)
        M0file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
    end
    
%% III-e. Lastly, still can't find it (error)
    
    if isempty(M0file)
        error('Cannot find the M0file')
    end
else
    M0file=[];
end

%
%%%
%%%%%
%%%
%

%% IV-a. Check for a brain mask

if whBM
    BMfile=[];
    % let's find the most up-to-date / full mask.
    if isfield(mrQ,'BrainMask')
        
        if  exist(mrQ.BrainMask,'file')
            BMfile= mrQ.BrainMask;
        end
    end
    if isfield(mrQ,'HMfile')
        
        if  exist(mrQ.HMfile,'file')
            BMfile= mrQ.HMfile;
        end
    end
    if isfield(mrQ,'FullMaskFile')
        
        if  exist(mrQ.FullMaskFile,'file')
            BMfile= mrQ.FullMaskFile;
        end
    end
    
    if isempty(BMfile)
        BMfile = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
    end
    
%% IV-b. Can't find it (error)
    
    if isempty(BMfile)
        error('Can not find the brain mask file')
    end
else
    BMfile=[];
end

%% Done