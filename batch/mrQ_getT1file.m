function T1file=mrQ_getT1file(mrQ)
% search for the T1 file. can be one of three file names depending on the
% fit process.


T1file=[];
%% is there an linear fit T1 file?
if isfield(mrQ.AnalysisInfo,'T1LFfile')

    if ( exist(mrQ.AnalysisInfo.T1LFfile,'file') == 0 )
        T1file= mrQ.AnalysisInfo.T1LFfile;
    end
end
%% is it a wieght linear t1 file
if isfield(mrQ.AnalysisInfo,'T1WLFfile')

    if exist(mrQ.AnalysisInfo.T1WLFfile,'file')
        T1file= mrQ.AnalysisInfo.T1WLFfile;
    end
end
%% is it an non linear t1 file?
if isfield(mrQ.AnalysisInfo,'T1lsqfile')

    if  exist(mrQ.AnalysisInfo.T1lsqfile,'file') 
        T1file= mrQ.AnalysisInfo.T1lsqfile;
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

%% select manually
   if isempty(T1file)
   T1file = mrvSelectFile('r','*.nii.gz','Select T1 fit file',mrQ.spgr_initDir);
   end
   
%% Last stll can't find it an error

if isempty(T1file)
    error('Can not find the T1file')
end

%% Done

   