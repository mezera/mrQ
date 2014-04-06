function [s coilNum] = makeStructFromNifti(niftiFile,multiChannels,struc,permution)
% 
% [s] = makeStracFromNifti(niftiFile,multiChannels,struc)
% 
% Load a nifti file and make a structure that the Align tool (knk's?) knows
% how to work with. The data needs to be reshaped and permuted so that data
% ---> (x,y,z,coils).
% 
% INPUTS: 
%       niftiFile     - full path to nifti file.
%       multiChannels - If multiChannels = -1 this will take all the
%                       channels and separate them. If it's a different
%                       number then we will work only on those channels
%                       named in multiChannels (vector). If it = -2 then we
%                       use only the last channel. *** EXPLAIN IN MORE
%                       DETAIL ***
%     struc             provide a structior that the data will be add to
%     permution           rearange the data. defult 1. most time the data from multi chanel need a permutation if this is not neccery use 0. 
% OUTPUTS:
%       [s]           - Structure containing the nifti data in a format
%                       that works well with the Align tool. 
% 
% SEE ALSO:
%   mrQ_initSPGR.m
% 
% 
% (C) Stanford University, VISTA Lab. 
% 

%% Check inputs and read in data

if ~exist('struc','var')
    struc = [];
end
if ~exist('permution','var') || isempty(permution)
    permution = 0;
end
DD = niftiRead(niftiFile);


%% Reshape and permute data - create the structure

if exist('multiChannels','var') && ~isempty(multiChannels)
    
    % Reshape the data structure so the 3rd and 4th dimensions are
    % switched-up and Permute the data so that dat dim --> (x,y,z,coils)
    sz  = size(DD.data);
    
    if length((sz))<4
        coilNum=1;
        sz(4)=1;
    else
    coilNum=sz(4)-1;
    end
    if permution==1
        dat = reshape(DD.data,[sz(1) sz(2) sz(4) sz(3)]);
        dat = permute(dat,[1 2 4 3]);
    else
        dat=DD.data;
    end
    % Use all coils (4D)
    if multiChannels == -1, 
        for i = 1:sz(4)
            s(i) = makeStruct(dat(:,:,:,i),DD,['coil ' num2str(i) ' from ' DD.fname],struc); %#ok<AGROW>
        end   
        
    % Use only the last coil (4D)    
    elseif multiChannels == -2 
        s(1) = makeStruct(dat(:,:,:,sz(4)),DD,['last coil from ' DD.fname],struc);

     % Use only the elected coils        
    else
        for i = 1:numel(multiChannels)
            s(i) = makeStruct(dat(:,:,:,multiChannels(i)),DD,['coil ' num2str(multiChannels(i)) ' from ' DD.fname],struc); %#ok<AGROW>
        end
    end
    
%for single chanel nifti
else
    s(1) = makeStruct(DD.data(:,:,:,end),DD,['image from ' DD.fname],struc);
coilNum=1;
end

return


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function s = makeStruct(dat,DD,seriesDesc,struc)
if (exist('struc','var') && ~isempty(struc))
    s(1)        = struc;
    s(1).imData = dat;
else
    s(1).imData            = dat;
    s(1).imToScanXform     = DD.sto_xyz;
    s(1).mmPerVox          = DD.pixdim(1:3);
    s(1).dims              = DD.dim(1:3);
    s(1).seriesDescription = seriesDesc;
end


