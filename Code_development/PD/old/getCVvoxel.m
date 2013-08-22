function [hold use] =getCVvoxel(Nvoxels,HoldForCV)
%[hold use] =getCVvoxel(Nvoxels,HoldForCV)
% this function randomly spleetset of voxel(location) form the pool of
% location (maximum location =Nvoxels) to two group. the hold and use
% group. the fraction on location in use defined by HoldForCV that need to
% be value betwen 0 to 1.

if notDefined('HoldForCV') 
    HoldForCV=0;
end

if HoldForCV<0 || HoldForCV>1
    HoldForCV=0;
end

% get a radom order of the voxels position
N=randperm(Nvoxels);

% define how many voxel to hold
Lastvoxel=round(Nvoxels*HoldForCV);

% the voxel position to hold
hold=N(1:Lastvoxel);

% the other voxel 
use=N(Lastvoxel+1:end);
