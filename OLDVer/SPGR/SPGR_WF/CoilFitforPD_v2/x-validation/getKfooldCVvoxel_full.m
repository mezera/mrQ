function [useX kFold] =getKfooldCVvoxel_full(Nvoxels,Ncoils,kFold)
%[hold use] =getKfooldCVvoxel_full(Nvoxels,Kfold)
% this function randomly spleetset of voxel(location) form the pool of
% location (maximum location =Nvoxels) to two group Kfold times. the hold and use
% group. the fraction of data n each group is defined by the Kfol0d 2 50%
% precent of the data. 10 Kfolod 0.1 10% ect.
% ....

if notDefined('kFold')
    kFold=1;
end

if kFold<0 || kFold>Ncoils
        fprintf(' kFold (= %d) is greater then the maxsimal posibale value. Changing to maxsimal kFold (=  %d) \n',kFold,Ncoils)
    kFold=Ncoils;
end

useX=ones(Nvoxels,Ncoils);
    % loop and get a random order of the coil to fit for each voxel

for ii=1:Nvoxels
    
    if  kFold==Ncoils;
        useX(ii,:)=randperm(kFold);
    end
    
    if  kFold<Ncoils;
        N=ceil(Ncoils/kFold);
        list=[];
        for jj=1:N;
            list=[list randperm(kFold)];
        end
        useX(ii,1:Ncoils)=list(1:Ncoils);
    end
    
    
end

