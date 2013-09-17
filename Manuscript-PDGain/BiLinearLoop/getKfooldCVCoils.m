function [useX ] =getKfooldCVCoils(Nvoxels,Ncoils,kFold) 



 useX=zeros(Nvoxels,Ncoils);
 % loop and get a random order of the coil to fit for each voxel
 
 for ii=1:Nvoxels
     useX(ii,:)= randperm(Ncoils);
 end
   

 
if  kFold<Ncoils;
        N=ceil(Ncoils/kFold);
        list=[];
        for jj=1:N;
            list=[list randperm(kFold)];
        end
        useX(ii,1:Ncoils)=list(1:Ncoils);
    end
    
    