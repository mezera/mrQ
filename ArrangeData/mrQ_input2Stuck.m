function [s, niifiles] = mrQ_input2Stuck(input,includeIm)
%
% [s, niifiles] = mrQ_input2Stuck(input,includeIm)
% 
%   This function takes the NIfTI files and the scan parameters, and puts
% them into a .mat structure which is used by the mrQ code. This is an
% alternative to reading the information from the DICOMs. Getting away from
% the DICOM is useful when dealing with different scanners' convention dcm
% files. Our code works smoothly only with GE DICOMs, so this is a way to
% solve the dcm problem.
%    Note that the user should be careful when inputting the scan
% parameters (TR, TE, flip angle, inversion time, field strength, and NIfTI
% name and directory).
% 
% INPUTS:
%   input      - A structure with the scanning parameters as provided
%                by the user. It may include these fields: dir (a directory
%                where the NIfTI is saved); name (NIfTI image name or part
%                of it); flip angle; TR; TE; IT (inversion time); or field
%                strength.
%   includeIm  - Determines whether to include the image from the NIfTI in
%                the output structure "s". (default is 1, "yes")
%
% OUTPUTS:
%   s          - A .mat strucute with the input information that can be
%                used by the mrQ code
%   niifiles   - The list of the NIfTI files that were used to make s.
%
% Example:
%   inputData_seir.dir  = mrQ.RawDir; 
%   inputData_seir.name = {'1301' '1401' '1501' '1601'}; 
%   inputData_seir.TR   = [3000 3000 3000 3000]; 
%   inputData_seir.TE   = [38 38 38 38]; 
%   inputData_seir.IT   = [2400 1200 400 50];
%   mrQ = mrQ_Set(mrQ,'inputdata_seir',inputData_seir)
%   s = mrQ_input2Stuck(mrQ.inputdata_seir);
% 
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
% 2015
%



if notDefined('includeIm')
    includeIm=1;
end


for i=1:length(input.name)
    
    % Note that input.name might already be the full path to the file. For
    % example, if you use mrQ_initInputData, they will all be full paths.
    if exist(input.name{i},'file')
        tmp=input.name{i};
    else
        cmd=fullfile(input.rawDir,[input.name{i} '*']); % BF input.dir = input.rawDir
        tmp=ls(cmd);
        tmp=tmp(1:end-1);
    end
    
    I = readFileNifti(tmp);
    niifiles{i} = tmp;
    
    if includeIm==1
        s(i).imData=I.data;
    end
    
    s(i).mmPerVox=I.pixdim;
    s(i).imToScanXform=I.qto_xyz;
    s(i).dims=I.dim;
    [~, name]=fileparts(I.fname);
    [~, name]=fileparts(name);
    
    s(i).seriesDescription=name;
end

% flip angle in degrees
if isfield(input,'flipAngle')
    for i=1:length(input.name)
        s(i).flipAngle=input.flipAngle(i);
    end;
end;

% TR in milliseconds
if isfield(input,'TR')
    for i=1:length(input.name)
        s(i).TR=input.TR(i);
    end;
end;

% TE in milliseconds
if isfield(input,'TE')
    for i=1:length(input.name)
        s(i).TE=input.TE(i);
    end;
end;

% IT in milliseconds
if isfield(input,'IT')
    for i=1:length(input.name)
        s(i).inversionTime=input.IT(i);
    end;
end;

%field strength
if isfield(input,'fieldStrength')
    for i=1:length(input.name)
        s(i).fieldStrength=input.fieldStrength(i);
    end;
end;

if isfield(input,'orientation')
    % if the convention of X, Y, and Z need to be flipped
    for ii=1:length(input.name)
        s(ii).imToScanXform(1,1)=input.orientation(1)*s(ii).imToScanXform(1,1);
        s(ii).imToScanXform(2,2)=input.orientation(2)*s(ii).imToScanXform(2,2);
        s(ii).imToScanXform(3,3)=input.orientation(3)*s(ii).imToScanXform(3,3);
        
    end
end

