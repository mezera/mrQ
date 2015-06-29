function [T1file, M0file,BMfile]=mrQ_get_T1M0_files(mrQ,whT1,whM0,whBM)
% [M0file, T1file]=mrQ_get_T1M0_files(mrQ)
%
% search for the T1 and M0 files. can be one of three file names depending on the
% fit process.

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
    
    
    %% is there an linear fit T1 file?
    %Find T1
    
    if isfield(mrQ,'T1_LFit')
        if ( exist(mrQ.T1_LFit,'file') == 0 )
            T1file= mrQ.T1_LFit;
        end
    end
    
    %  the same but with larger  mask
    if isfield(mrQ,'T1_LFit_HM')
        if ( exist(mrQ.T1_LFit_HM,'file') == 0 )
            T1file= mrQ.T1_LFit_HM;
        end
    end
    
    %  Maybe after B1 corraction
    if isfield(mrQ,'T1_B1_LFit')
        if ( exist(mrQ.T1_B1_LFit,'file') == 0 )
            T1file= mrQ.T1_B1_LFit;
        end
    end
    
    %% is it a wieght linear t1 file
    if isfield(mrQ,'T1WLFfile')
        
        if exist(mrQ.T1WLFfile,'file')
            T1file= mrQ.T1WLFfile;
        end
    end
    
    %% is it an non linear t1 file?
    if isfield(mrQ,'T1_B1_lsqFit')
        
        if  exist(mrQ.T1_B1_lsqFit,'file')
            T1file= mrQ.T1_B1_lsqFit;
        end
    end
    
    %% select manually
    if isempty(T1file)
        T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
    end
    
    %% Last still can't find it an error
    
    if isempty(T1file)
        error('Can not find the T1file')
    end
end
%% is there an linear fit M0 file?
%Find M0
if whM0
    M0file=[];
    if isfield(mrQ,'M0_LFit')
        if ( exist(mrQ.M0_LFit,'file') == 0 )
            M0file= mrQ.M0_LFit;
        end
    end
    
    %  the same but with larger  mask
    if isfield(mrQ,'M0_LFit_HM')
        if ( exist(mrQ.M0_LFit_HM,'file') == 0 )
            M0file= mrQ.M0_LFit_HM;
        end
    end
    
    
    %  Maybe after B1 corraction
    if isfield(mrQ,'M0_B1_LFit')
        if ( exist(mrQ.M0_B1_LFit,'file') == 0 )
            M0file= mrQ.M0_B1_LFit;
        end
    end
    
    
    %% is it a wieght linear M0 file
    
    if isfield(mrQ,'M0_B1_LWFit')
        
        if exist(mrQ.M0WLFfile,'file')
            M0file= mrQ.M0WLFfile;
        end
    end
    %% is it an non linear M0 file?
    
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

    %% maybe it is already in the outputDir
    if isfield(mrQ, 'maps')
        if isfield(mrQ.maps,'T1path')
            if exist(mrQ.maps.T1path,'file')
                T1file= mrQ.maps.T1path;
            end
        end
    end
    
    % Note that we are not move M0 image to the output Dir so we are not looking for it there.
    
    
    
    if isempty(M0file)
        M0file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
    end
    
    %% Last still can't find it an error
    
    if isempty(M0file)
        error('Can not find the M0file')
    end
end
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
    
    %% Last still can't find it an error
    
    if isempty(BMfile)
        error('Can not find the M0file')
    end
end

%% Done

