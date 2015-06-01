function [s, niifiles] = mrQ_input2Stuck(input,includeIm)
%
% [s, niifiles] = mrQ_input2Stuck(input,includeIm)
% 
% This function take the nifti files and the scan parameters and put them
% in mat structure that is used by the mrQ code. this is an alternative to
% reading the information from the dicom. getting away from the dicom is
% usful when dealing with different scanner convention dcm file. our code
% is working smoothly only with GE dicom so this is a way to solve the dcm
% problem. not that the user have to be carfule when plaging the scan
% parameters (TR ,TE flipangle, inversion time field Strength and nifti
% name and dir)
% 
% INPUTS:
%   input -     a structure with the scanning parameters that is provide
%               by the user and may include these feileds [dir (a directory were
%               the niffti saved); name (Niifti image name or part of it);
%               flipAngle; TR; TE; IT; (inversionTime); fieldStrength]
%   includeIm-  to include the image from the nifti in the structure (s)
%   includeIm   = 1 (defult).
%
% OUTPUTS:
%   s -         A mat strucute with the input information that can be used
%               by mrQ code 
%   niifiles -  The list of the nifti files that was used to make s
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


if notDefined('includeIm')
    includeIm=1;
end


for i=1:length(input.name)
    % input.name might already the full path to the file. For example, if
    % you use mrQ_initInputData they will all be full paths.
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

% flipAngle in degree
if isfield(input,'flipAngle')
    for i=1:length(input.name)
        s(i).flipAngle=input.flipAngle(i);
    end;
end;

% TR in milisec
if isfield(input,'TR')
    for i=1:length(input.name)
        s(i).TR=input.TR(i);
    end;
end;

% TE in milisec
if isfield(input,'TE')
    for i=1:length(input.name)
        s(i).TE=input.TE(i);
    end;
end;

% IT in milisec
if isfield(input,'IT')
    for i=1:length(input.name)
        s(i).inversionTime=input.IT(i);
    end;
end;
if isfield(input,'fieldStrength')
    for i=1:length(input.name)
        s(i).fieldStrength=input.fieldStrength(i);
    end;
end;

if isfield(input,'orientation')
    % if the convention of X Y and Z need to be fliped
    for ii=1:length(input.name)
        s(ii).imToScanXform(1,1)=input.orientation(1)*s(ii).imToScanXform(1,1);
        s(ii).imToScanXform(2,2)=input.orientation(2)*s(ii).imToScanXform(2,2);
        s(ii).imToScanXform(3,3)=input.orientation(3)*s(ii).imToScanXform(3,3);
        
    end
end

