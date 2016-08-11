function [d, xform,AlignDir ]=mrQ_fslAlignCall(Dpath,d,niilist,Ref)
%[d, xform,AlignDir]=mrQ_fslAlignCall(Dpath,d,niilist,Ref)
%
%  This function is a wrapper that registers a niifile in a directory 
% 'Dpath' to a reference file. It uses the fsl command 'flirt', and the 
% call uses rigid body. The output "d" is a matlab structure which contains 
% all of the the images and meta information from the nifti files. 
% If "d" is also input, it will add to the existing structure. 
% AlignDir is where the alignment nifti files are saved.
%
% INPUTS:
% >> Dpath is the directory of the reference files; its default is pwd.
% >> d is the structure which will contain all the images and meta
% information.
% >> niilist is the list of nifti files.
% >> Ref is the filename for the reference files; its default is the first
% nifti from the niilist.
%
% OUTPUTS:
% >> d, the new/changed structure which will contain all the images and meta
% information.
% >> xform, the transformation matrix for the alignment.
% >> AlignDir, the directory where the aligned nifti files are saved.
%
% AM vista lab 2013
%(C) Stanford University, VISTA
% Edited by Shai Berman and Jonathan Bain, June-02-2015
%

if notDefined('Dpath')
    Dpath=pwd;
end

if notDefined('niilist')
niilist=dir(fullfile(Dpath,'*nii.gz'));
end

if notDefined('Ref') 
    
    % Use the first nifti to register, in case Ref is not defined
%     Ref=fullfile(Dpath,niilist(1).name); >>>>shai
Ref=fullfile(niilist{1});
end

if notDefined('d') 
d=[];
end

AlignDir=fullfile(Dpath,'Align');

if  ~exist(AlignDir,'dir');mkdir (AlignDir);end

for ii=1:length(niilist)
    
%     inPutFile=fullfile(Dpath,niilist(ii).name); >>>>>shai
     inPutFile=fullfile(niilist{ii}); 
        outPutFile=fullfile(AlignDir,['Align_ImNo' num2str(ii) '.nii.gz']);
        outPutmat=fullfile(AlignDir,['Align_ImNo' num2str(ii) '.mat']);

    eval(['! flirt -in ' inPutFile ' -ref ' Ref  ' -out ' outPutFile  ' -omat ' outPutmat ' -dof 6'])
    
    Alignimage=readFileNifti(outPutFile);
    d(ii).imData=Alignimage.data;
        d(ii).imToScanXform=Alignimage.qto_xyz;
        % NOTE: the xform comes in meters. Multiply by 1000 to get millimeters. 
   d(ii).imToScanXform(1:3,1:3)=d(ii).imToScanXform(1:3,1:3)*1000;
    d(ii).mmPerVox=Alignimage.pixdim*1000;
    
    
end

 xform=   d(1).imToScanXform;