function [holdX useX] =getKfooldCVvoxel(Nvoxels,Kfolod)
%[hold use] =getKfooldCVvoxel(Nvoxels,Kfold)
% this function randomly spleetset of voxel(location) form the pool of
% location (maximum location =Nvoxels) to two group Kfold times. the hold and use
% group. the fraction of data n each group is defined by the Kfol0d 2 50%
% precent of the data. 10 Kfolod 0.1 10% ect.
% ....

if notDefined('Kfolod') 
    Kfolod=1;
end

if Kfolod<0 || Kfolod>Nvoxels
    Kfolod=1;
end

% get a radom order of the voxels position
N=randperm(Nvoxels);

locs= round(  linspace(1,Nvoxels,Kfolod+1)); 

useX=ones(Nvoxels,Kfolod);
holdX=zeros(Nvoxels,Kfolod);

for ii=1:Kfolod
    
holdX(N(locs(ii):locs(ii+1)),ii)=1;
useX(N(locs(ii):locs(ii+1)),ii)=0;
end