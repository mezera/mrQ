function [useX ] =getSparceCVvoxel_full(Nvoxels,Ncoils)
%
 %[useX ] =getSparceCVvoxel_full(Nvoxels,Ncoils)



 useX=zeros(Nvoxels,Ncoils);
 % loop and get a random order of the coil to fit for each voxel
 
 for ii=1:Nvoxels
     N=    randperm(Ncoils);
     useX(ii,N(1))=1;
 end
   

