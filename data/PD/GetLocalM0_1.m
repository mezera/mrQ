function [M01 SZ meanVal XX YY ZZ]=GetLocalM01(M0,boxNeibhors,Cvoxel)
%
%
% Obsolete - delete me
%
%

XX(1)=Cvoxel(1)+boxNeibhors(1,1);
XX(2)=Cvoxel(1)+boxNeibhors(2,1);
YY(1)=Cvoxel(2)+boxNeibhors(1,2);
YY(2)=Cvoxel(2)+boxNeibhors(2,2);
ZZ(1)=Cvoxel(3)+boxNeibhors(1,3);
ZZ(2)=Cvoxel(3)+boxNeibhors(2,3);
M01=M0(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2),:);
    
SZ=size(M01);
   [meanVal Incoil]=sort(squeeze(mean(mean(mean(M01)))),'descend');
M01=M01(:,:,:,Incoil);

   [meanVal Incoil]=sort(squeeze(mean(mean(mean(M01)))),'descend');

end
