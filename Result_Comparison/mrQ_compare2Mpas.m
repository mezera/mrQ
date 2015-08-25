function mrQ_compare2Mpas(T1wfile1_LR,T1wfile2_HR,outPutDir,morefiles2_HR,InsavefileN)
%
%function mrQ_compare2Mpas(file1,file2,otherFiles1,otherFile2)
% This function will register two T1 weighted files (or any other two) of the
% same subjects (coming from different scan or scanner) and aligns them
% using semi-automatic tools (devloped by Kendrick Kay). We use only rigid
% or affine registration parameters.


keyboard

% load the two volume
im1=readFileNifti(T1wfile1_LR);
im2=readFileNifti(T1wfile2_HR);

% define the region to align on. manually
[f,mn,sd] = defineellipse3d(double(im1.data));


% look and start the alinment manualy
alignvolumedata(im2.data,im2.pixdim,im1.data,im1.pixdim);

% to re-define:
%   [f,mn,sd] = defineellipse3d(double(im1.data),[],[],mn,sd);
%   alignvolumedata(im2.data,im2.pixdim,im1.data,im1.pixdim,ttr);


for i=1:10
    i
 alignvolumedata_auto(mn,sd,0,[1 1 1]); 
  alignvolumedata_auto(mn,sd,1,[1 1 1]);

alignvolumedata_auto(mn,sd,2,[1 1 1]);
end;

ttr=alignvolumedata_exporttransformation;
warpParmFile= fullfile(outPutDir,'T1w2_to_T1w1ManParam');
save(warpParmFile,'ttr','md','sd')

close all

if notDefined('InsavefileN');InsavefileN=[];end
if notDefined('morefiles2_HR');morefiles2_HR=[];end

mrQ_knkREgisterIm(T1wfile1_LR,T1wfile2_HR,warpParmFile,outPutDir,morefiles2_HR,InsavefileN)
    