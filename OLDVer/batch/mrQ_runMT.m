function mrQ_runMT(mrQ)
%RawDir,arrangeRawFlag,MT_seriesNumbers,MT_offsets,sub,refIm,interp,mmPerVox,skip,clobber)

%mrQ_runMT(mrQ,RawDir,arrangeRawFlag,MT_seriesNumbers,refIm,sub,clobber)
%
% # a batch call to run all the mrQ MT fits in one click
% # the batch  run over mrQ functions to claculate the maps MT's (BPF,K) it will also run SGE and softwhere. this and
%%the mrVista reposotorry must be set and updated

%   RawDir            - The path of the directory containing the raw MR
%                       data.the function assume the data are arrange in the same way as it is in the Stanford CNI magnet.
%                       the data should be saved as niifti files and a dicom directory with dicom are also exist.
%                       This convention may be different for each scanner.
%                       It is good idea to cheack this and maybe arrange the data with this conventions or edit the first function that call the raw directory ( mrQ_arrangeData.m)
%
%   arangeRawFlag     - Needs to be 1 to move every file from the original
%                       magnet directory (dataDir) to a raw directory
%                       inside it (if this was done set = 0, or empty, defult).
%
%   MT_seriesNumbers- A cell array containing the MT series numbers
%                       from dataDir/raw.  (if this input exists a MT dir
%                       will be created with data  directories
%                       inside it). The a smple of the scan dicom will be extracted from
%                       the data Dir as well as the nifti. If MT dit already exists and a new
%                       one will be created.
%
%   MT_seriesNumbers- A vector with the MT off sets.
%                       this input need to be given as the offset are not save in the dicoms or nifti.
%                       tipacly we save this information in the name of the scan
%
%  refImg              Different ref images can be used as an input (refImg is a
%                      path to a nifti image).
%
%       sub            - Subject name used for different file names
%                       like the SGE. if the name is not porvided we will take it for the
%                        raw directory
%
%  clobber              -Overwrite existing data and reprocess. [default = false]

%output % save a mrQparm file the document the anaysis. the different file
%make fiels and directory. the run end when a map directory with the map
%is crated including f and k maps files
%see also mrQ_run
%Example :
%
%
%
%

%%

if ~isfield(mrQ,'rawDir')
    disp('Raw data need to be provided');
    error
end;

% % do we like to arange the data or it was done before
if ~isfield(mrQ,'arrangeRawFlag')
    mrQ.arrangeRawFlag=0;
    % this will be zerow if a raw dir was allready made and data was moved
    % before
end;


%name=fullfile(mrQ.RawDir,'mrQ_params.mat');

% if (exist(name,'file'))
%     load(name);
%      fprintf('load mrQ data !              \n ');
%     
% end

if isfield(mrQ,'arangingMT_done');
else
    mrQ.arangingMT_done=0;
end
    if (mrQ.arangingMT_done==0);


    fprintf('aranging the MT data \n ')
     [mrQ.MT]=  mrQ_arrangeDataMT(mrQ.rawDir,mrQ.arrangeRawFlag,mrQ.MT.MT_seriesNumbers)
        fprintf('arange data is done - done! \n');
     mrQ.arangingMT_done=1;
    save(mrQ.name,'mrQ');
    else
             fprintf(['using  data in ' mrQ.MT.MTDir               '\n ']);
    end
%%
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load(name);

if isfield(mrQ.MT,'MT_init');
else
    mrQ.MT.MT_init=0;
end

if     mrQ.MT.MT_init==0
    if ~isfield(mrQ.MT,'interp')
        
        mrQ.MT.interp=mrQ.interp;
    end
    if ~isfield(mrQ.MT,'mmPerVox')
       mrQ.MT.mmPerVox =mrQ.mmPerVox;
    end;
    if ~isfield(mrQ.MT,'refIm')
        mrQ.MTrefIm=mrQ.refIm;
    end
    if ~isfield(mrQ.MT.skip,'refIm')
        mrQ.MT.skip=mrQ.skip;
    end;
    %keep track of the variable we use for detail see inside the function
    
    [~, ~, ~,~,mrQ]=mrQ_initMT(mrQ.MT.MTDir,mrQ.MT.refIm,mrQ.MT.mmPerVox,mrQ.MT.interp,mrQ.MT.skip,mrQ.MT.MT_offsets,[],mrQ);
    mrQ.MT.MT_init=1;
    
    save(mrQ.name,'mrQ');
    fprintf('init MT - done!              \n ');
else
    fprintf('load init MT data              \n  ');
    
end
    
%    load(name);

if isfield(mrQ.MT,'MT_T1fit');
else
    mrQ.MT.MT_T1fit=0;
end

if     mrQ.MT.MT_T1fit==0
    
    mrQ_FIT_MT(mrQ)
    
end
    
    
    
    