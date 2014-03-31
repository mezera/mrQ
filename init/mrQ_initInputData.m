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
fprintf('[%s]\n',mfilename);

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
disp('Reading parameters from nifti files');
for jj = 1:numel(niFiles)
    nifti{end+1} = niftiGetParamsFromDescrip(niFiles{jj});
end

% Remove empty entries 
nifti(cellfun(@isempty,nifti)) = []; 


%% Flip Angle

% Rules: 
%	- Data must be 4D of same size (multi-coil) - don't know about the same size thing? 
%	- Minimum of 2 unique flip angles with 1 <=10 and 1 >= 10. 
%	- FA must be <=30
%	- TE <= 3
% 	- TI == 0 
%	- OFF RESONANCE == 0

disp('Looking for nifti files with valid flip angles...');
fa = {};
for ii = 1:numel(nifti)
	if isfield(nifti{ii},'fa') && ( isfield(nifti{ii},'rs') || isfield(nifti{ii},'r') )
		if ( nifti{ii}.fa <= 30 ) && ( nifti{ii}.te <= 3 ) && ( nifti{ii}.ti == 0 )
			fa{end+1} = ii;
		end
	end
end

fprintf('%d nifti files found:\n',numel(fa));
for nn = 1:numel(fa)
	[p, f, e ] = fileparts(nifti{fa{nn}}.niftiFile); 
	thename = [f, e];
	epoch = explode('/',p); 
	fprintf('\t%s/%s: \tFlip angle = %d\n',epoch{end},[f,e],nifti{fa{nn}}.fa);
end

% Is 1 flip angle less than 10 and one greater than 10
f = zeros(size(fa));
for fas = 1:numel(fa)
	f(fas) = nifti{fa{fas}}.fa;
end
if numel(f) < 3 && ( isempty(find(f > 10)) || isempty(find(f < 10)) ) 
	warning('There are fewer than 3 valid nifti files and those flip angles are not consistent with this analysis!');
end
			

%% Inversion Time
% RULES: 
%	- 3 unique TIs 
%	- TI must be non-zero
% 	- TI must be <=2500 * 
%	- Constant TE

disp('Looking for nifti files with valid inversion times...');
it ={};
for jj = 1:numel(nifti)
    if isfield(nifti{jj},'ti')
    	if ( nifti{jj}.ti ~=0 ) && (nifti{jj}.ti <= 2500)
    		it{end+1} = jj;
        end
    end
end

fprintf('%d nifti files found:\n',numel(it));
for nn = 1:numel(it)
	[p, f, e ] = fileparts(nifti{it{nn}}.niftiFile); 
	epoch = explode('/',p); 
	fprintf('\t%s/%s: \tInversion time = %d\n',epoch{end},[f,e],nifti{it{nn}}.ti);
end

% Find out if the inversion times are unique
t = zeros(size(it));
for nit = 1:numel(it)
	t(nit) = nifti{it{nit}}.ti;
end
if numel(unique(t)) ~= numel(it)
	warning('INVERSION TIMES ARE NOT UNIQUE!');
end

% Find out if the TEs are constant
e = zeros(size(it));
for ne = 1:numel(it)
	e(ne) = nifti{it{ne}}.te;
end
if numel(unique(e)) ~= 1
	warning('ECHO TIMES ARE NOT CONTSTANT!');
end


%% Construct inputData structure
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

if numel(it) < 3 || numel(fa) < 2
    fprintf('[%s] - Files required for mrQ not found in %s\n',mfilename,mrQ.RawDir);
    else
    fprintf('[%s] - Files required for mrQ found!\n',mfilename);
end

return
