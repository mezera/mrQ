function PDmedianfromGain(M0ftfull,xform,outDir,j)
 V=logical(M0ftfull); 
 V1=sum(V,4);
  V1=logical(V1);
  tmp=zeros(size(V1));
  wh=find(V1);
  for q=1:length(wh)
    [xx yy zz]=  ind2sub(size(V1),wh(q));
   tmp(xx,yy,zz) =   median(nonzeros(M0ftfull(xx,yy,zz,:)));
   q
  end;
    WFfile=fullfile(outDir,['PDsqrt_fitboxMedian_' num2str(j) '.nii.gz']);

 dtiWriteNiftiWrapper(single(tmp), xform, WFfile);
   %             keyboard;