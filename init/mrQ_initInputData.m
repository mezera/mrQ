function mrQ = mrQ_initInputData(mrQ)
% 
%  mrQ = mrQ_initInputData(mrQ)
% 
% This function will take a mrQ structure with a 'RawDir' field and
% recursively parse all nifti files in that directory, returning a mrQ
% structure with inputdata_spgr and inputdata_seir structures attached.
% 
%* This function only works with data that were run through NIMS.
% 
% Roughly how it works:
%   1. Finding all the nifti files in mrQ.RawDir
%   2. Extracting the params for each nifi file (using
%      niftiGetParamsFromDescrip.m)
%   3. Sorting the nifti files based on their params into two structures
%      (1 for inputdata_spgr and 1 for inputdata_seir).
%   4. It uses the Flip Angle and the Inversion Times of each nifti file to
%      determine in which structure each nifit file belongs. 
% 
% EXAMPLE USAGE:
%   mrQ = mrQ_Create('/biac4/wandell/data/westonhavens/upload/testlab/20130509_1152_4534')
%   mrQ = mrQ_initInputData(mrQ);
%   mrQ = 
%                      RawDir: [1x66 char]
%                      ...
%              inputdata_spgr: [1x1 struct]
%              inputdata_seir: [1x1 struct]
% 
%   mrQ.inputdata_seir =
% 
%           name: {[1x109 char]  [1x110 char]  [1x110 char]  [1x108 char]}
%             TR: [1 1 1 1]
%             TE: [47 47 47 47]
%             IT: [400 2400 1200 50]
%         rawDir: '/biac4/wandell/data/westonhavens/upload/testlab/20130509_1152_4534'
% 
%   mrQ.inputdata_spgr =
% 
%         name: {[1x108 char]  [1x106 char]  [1x108 char]  [1x107 char]}
%                TR: [33 33 33 33]
%                TE: [2 2 2 2]
%         flipAngle: [20 10 30 4]
%            rawDir: '/biac4/wandell/data/westonhavens/upload/testlab/20130509_1152_4534'
% 
% 
% 
% (C) VISTA Lab, Stanford University 2014
% 


%% Process inputs and checks
% mrQ.RawDir = '/biac4/wandell/data/westonhavens/upload/testlab/20130509_1152_4534';

if ~isfield(mrQ,'RawDir')
    error('No rawDir in mrQ structure');
end

if ~isfield(mrQ,'fieldstrength')
    mrQ = mrQ_Set(mrQ,'fieldstrength',3);
end



%% Build a struct with the path to each nifti file in mrQ.RawDir

% Find all nifti files in mrQ.RawDir (recursively)
tn = tempname;
cmd = ['find ' mrQ.RawDir ' -follow -type f -name "*.nii.gz" | tee ' tn];

% Run the command 
[status, result] = system(cmd);
if status ~= 0 
    error('There was a problem finding nifti files.'); 
end

% niFiles will now have a full-path list of all nifti files
if ~isempty(result)
    niFiles = readFileList(tn);
else
    fprintf('mrQ - No nifti files found in %s', mrQ.RawDir); 
%     Should this be returned empty, or should some flag be set??
%     mrQ = []; 
    return
end


%% Read each nifti and set up the structures

nifti = {};

% Read the parameters for each of the nifti files
for jj = 1:numel(niFiles)
%     fprintf('Reading  %s\n',niFiles{jj});
%     if jj == 9
%         keyboard
%     end
    nifti{end+1} = niftiGetParamsFromDescrip(niFiles{jj});
end

% Remove empty entries 
nifti(cellfun(@isempty,nifti)) = []; 

% Get the flip angles 
fa = {};
for ii = 1:numel(nifti)
    if isfield(nifti{ii},'fa') && ( isfield(nifti{ii},'rs') || isfield(nifti{ii},'r') )
        switch nifti{ii}.fa
            case {10, 20, 30, 4};
                fa{end+1} = ii;
        end
    end
end

% Do another check of those nifti files for those with duplicate FAs


% Get the IT for the SEIR structure
it ={};
for jj = 1:numel(nifti)
    if isfield(nifti{jj},'ti')
        switch nifti{jj}.ti
            case {2400, 1200, 400, 50}
              it{end+1} = jj;
        end
    end
end

% Do another check of those nifti files for those with duplicate ITs



%% Construct inputData structures

inputData_spgr = {};

% Loop over the fa and it structs and populate the input data fields
for x = 1:numel(fa)
    inputData_spgr.name{x}          = nifti{fa{x}}.niftiFile;
    inputData_spgr.TR(x)            = nifti{fa{x}}.tr;
    inputData_spgr.TE(x)            = nifti{fa{x}}.te;
    inputData_spgr.flipAngle(x)     = nifti{fa{x}}.fa;
    inputData_spgr.fieldStrength(x) = mrQ.fieldstrength;
end

inputData_seir = {};
for y = 1:numel(it)
    inputData_seir.name{y}     = nifti{it{y}}.niftiFile;
    inputData_seir.TR(y)       = nifti{it{y}}.tr;
    inputData_seir.TE(y)       = nifti{it{y}}.te;
    inputData_seir.IT(y)       = nifti{it{y}}.ti;
    inputData_seir.orientation = [1 1 1]; % how can this be derived???
end


%% Return the results

inputData_seir.rawDir = mrQ.RawDir;
inputData_spgr.rawDir = mrQ.RawDir; 

mrQ = mrQ_Set(mrQ,'inputdata_spgr',inputData_spgr);
mrQ = mrQ_Set(mrQ,'inputdata_seir',inputData_seir);

if numel(it) <2 || numel(fa) <2
    fprintf('\n[%s] - Files required for mrQ not found in %s\n',mfilename,mrQ.RawDir);
end

return
