function mrQ_makeNiftiFromStruct(s,saveFileName,xform)
% 
% makeNiftiFromStruct(s,saveFileName,[xform])
% 
% Take a structure and make a nifti that the Align tool knows to work with.
% this is funny but needed. 
% 
% INPUTS:
%       s:            structure
%       saveFileName: name for the resulting nifti - should be a full path.
%       xform:        transform. Defaults to s(1).imToScanXform.
% 
% 
% (C) Stanford University, VISTA Lab
% 
% 

%%

if (~exist('xform','var')|| isempty(xform)),
    xform=s(1).imToScanXform;
end

if size(s)>1,    
    Dat = zeros([s.dim size(s)]);
    
    for i=1:size(s)
        Dat(:,:,:,i)=s(i).imData;
    end
    
    s(1).imData = Dat;  
end

 dtiWriteNiftiWrapper(single(s(1).imData),xform, saveFileName);
 
 return