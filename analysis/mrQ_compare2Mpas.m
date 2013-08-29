function mrQ_compare2Mpas(T1wfile1_LR,T1wfile2_HR,outPutDir,otherFile2)
%
%function mrQ_compare2Mpas(file1,file2,otherFiles1,otherFile2)
% this function will register two T1 wighter file (or any other two) of the
% same subjects comming form different scan or scanner) and aligh them
% using semi automatic tools (devloped by Kendric Kay). we use only rigid
% or affine registration parameters.


keyboard

% load the two voulume
im1=readFileNifti(T1wfile1_LR);
im2=readFileNifti(T1wfile2_HR);

% define the region to align on. manualy
[f,mn,sd] = defineellipse3d(double(im1.data));


% look and start the alinment manualy
alignvolumedata(im2.data,im2.pixdim,im1.data,im1.pixdim);

% to re define:
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
warpFile= fullfile(outPutDir,'WarpMan_T1w2_to_T1w1.nii.gz');

p=pwd; cd '~avivm/matlab/vistasoft/trunk/kendrick/kendrick/alignvolumedata/private'


ok=reslicevolume(0,ttr,'cubic',3,[],1,0,double(im2.data),double(im2.pixdim),size(im2.data),double(im1.data),double(im1.pixdim),[size(im1.data) 1]);

dtiWriteNiftiWrapper(single(ok),im1.qto_xyz,warpFile);

for d=1:length(morefiles2_HR)

    im2=readFileNifti(morefiles2_HR{d});
     file=dir(morefiles2_HR{d});
savefileN=fullfile(outPutDir,['WarpMan_' file.name]);
    ok=reslicevolume(0,ttr,'cubic',3,[],1,0,double(im2.data),double(im2.pixdim),size(im2.data),double(im1.data),double(im1.pixdim),[size(im1.data) 1]);
dtiWriteNiftiWrapper(single(ok),im1.qto_xyz,savefileN);
clear ok 

end

cd (p)

    