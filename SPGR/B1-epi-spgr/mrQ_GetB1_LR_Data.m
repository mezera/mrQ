function [S, t1, BM1,boxSize,UseVoxN, skip ,f1, XX, YY, ZZ]= mrQ_GetB1_LR_Data(opt,Res,BM,ind)
% Load  data from the M0 and T1 file for a local regreation cdefined by voxel index and opt
%
%    [S, t1, BM1,boxSize, skip ,f1, XX, YY, ZZ]= mrQ_GetB1_LR_Data(opt,Res,BM,ind)
%
%
% Helper routine to get the data


SZ=size(BM);
S=[]; t1=[]; BM1=[];boxSize=[]; skip=0;f1=[]; XX=[];YY=[]; ZZ=[];

[xi, yi, zi]=ind2sub(SZ,ind);
Bzide= round(opt.FilterSize./opt.pixdim);

if any(Bzide)==0, Bzide(find(Bzide==0))=1;end % no zerow box size
     
 XX(1)=xi-Bzide(1); XX(2)=xi+Bzide(1);
 YY(1)=yi-Bzide(2); YY(2)=yi+Bzide(2);
 ZZ(1)=zi-Bzide(3); ZZ(2)=zi+Bzide(3);



if XX(1)<1; XX(1)=1;  end;  if XX(2)>SZ(1); XX(2)=SZ(1);end
if YY(1)<1; YY(1)=1;  end;  if YY(2)>SZ(2); YY(2)=SZ(2);end
if ZZ(1)<1; ZZ(1)=1;  end;  if ZZ(2)>SZ(3); ZZ(2)=SZ(3);end



 if ~isempty(Res)   
    %get the location of the box we work on in image space (x,y,z)

    % Pull out the data
% SEIR T1
    t1=Res{1}.im(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2));

%Brain mask
BM1=logical(BM(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2)));
   %% if t1 values are worng we will skip the voxels
    Bad=isnan(t1) | isinf(t1) | t1==0;
        BM1(Bad)=0;
%raw SPGR
for ii=3:length(Res)
 tmp= Res{ii}.im(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2));
 Bad=isnan(tmp) | isinf(tmp) | tmp==0;
        BM1(Bad)=0;
S(:,:,:,ii-2)=tmp;
end
% This is the 4D size  of the box
boxSize = size(S);

   [f1] = makegaussian3d(boxSize(1:3),[0.5 0.5 0.5],[0.25 0.25 0.25]);
UseVoxN= length(find(BM1));
if UseVoxN<length(BM1(:))*opt.pracent_coverage ||  UseVoxN<10  % not enghf voxels
    skip=1;
end
 end



 
 

