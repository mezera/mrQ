function [d, xform,AlignDir ]=mrQ_fslAlignCall(Dpath,d,niilist,Ref)
%[d, xform,AlignDir ]=mrQ_fslAlignCall(Dpath,d,Ref)
%
%  a raper that register a niifile in a directory Dpath to a reference
%  file. using fsl flirt. the call use rigid body.
% the out put stracture d. is a matlab stracture that contain all the 
% the images and meta information from the niffti. if d is also input it
% will add on an exsisting stracture.
%AlignDir is where the alignment nifti are saved.
%
% AM vista lab 2013
%(C) Stanford University, VISTA
%
%

if notDefined('Dpath')
    Dpath=pwd;
end

if notDefined('niilist')
niilist=dir(fullfile(Dpath,'*nii.gz'));
end

if notDefined('Ref') 
    % use the first niifti to regiser in case res is not defined
    Ref=fullfile(Dpath,niilist(1).name);
    st=2;
else
    st=1;
end
if notDefined('d') 
d=[];
end

    AlignDir=fullfile(Dpath,'Align');

if  ~exist(AlignDir,'dir');mkdir (AlignDir);end

for ii=st:length(niilist)
    
    inPutFile=fullfile(Dpath,niilist(ii).name);
        outPutFile=fullfile(AlignDir,['Align_ImNo' num2str(ii) '.nii.gz']);
        outPutmat=fullfile(AlignDir,['Align_ImNo' num2str(ii) '.mat']);

    eval(['! flirt -in ' inPutFile ' -ref ' Ref  ' -out ' outPutFile  ' -omat ' outPutmat ' -dof 6'])
    
    Alignimage=readFileNifti(outPutFile);
    d(ii).imData=Alignimage.data;
        d(ii).imToScanXform=Alignimage.qto_xyz;
        % NOTE: the xfrom come in meter. i multipal by 1000 for mm. i hope
        % this is right. This is just a hake. this should be monitor 
   d(ii).imToScanXform(1:3,1:3)=d(ii).imToScanXform(1:3,1:3)*1000;
    d(ii).mmPerVox=Alignimage.pixdim*1000;
    
    
end
 xform=   d(1).imToScanXform;