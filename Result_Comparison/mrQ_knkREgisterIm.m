function mrQ_knkREgisterIm(T1wfile1_LR,T1wfile2_HR,KNK_ManParam,outPutDir,morefiles2_HR,InsavefileN)
%mrQ_knkREgisterIm(T1wfile1_LR,outPutDir,otherFile2)
% this function work with Kendrick tool for registration see knk git and
% svn reposatory. the function will register the nii file accordig to
% parameters save 

% load the two voulume
im1=readFileNifti(T1wfile1_LR);
im2=readFileNifti(T1wfile2_HR);
% and ther registration params
load(KNK_ManParam);
p=pwd; cd '~avivm/matlab/vistasoft/trunk/kendrick/kendrick/alignvolumedata/private'

if notDefined('InsavefileN')
    
    savefileN= fullfile(outPutDir,'WarpMan_T1w2_to_T1w1.nii.gz');
else
    savefileN=  InsavefileN{1};
end

ok=reslicevolume(0,ttr,'cubic',3,[],1,0,double(im2.data),double(im2.pixdim),size(im2.data),double(im1.data),double(im1.pixdim),[size(im1.data) 1]);

dtiWriteNiftiWrapper(single(ok),im1.qto_xyz,savefileN);
%register any other images (that have similar cordinate as T1wfile2_HR
for d=1:length(morefiles2_HR)

    im2=readFileNifti(morefiles2_HR{d});
     file=dir(morefiles2_HR{d});
     if notDefined('InsavefileN')
        
savefileN=fullfile(outPutDir,['WarpMan_' file.name]);
     else
       savefileN=  InsavefileN{d+1};
     end
    ok=reslicevolume(0,ttr,'cubic',3,[],1,0,double(im2.data),double(im2.pixdim),size(im2.data),double(im1.data),double(im1.pixdim),[size(im1.data) 1]);
dtiWriteNiftiWrapper(single(ok),im1.qto_xyz,savefileN);
clear ok 

end

cd (p)