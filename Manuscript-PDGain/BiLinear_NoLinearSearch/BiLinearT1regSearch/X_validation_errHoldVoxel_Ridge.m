function [err_X] = X_validation_errHoldVoxel_Ridge(g,M0,pBasis,nPositions,nCoils,Xmask)

%estimate coil coefficients across the volume
G = pBasis*g;
err_X=[];
% those are the location we need to Xvalidate with
loc=find(Xmask==1); 

nVoxels=length(loc);
M0_XV=M0(loc,:);
% we will spleat the X validation position to two. half for PD estimate and
[useX kFold] =getKfooldCVvoxel_full(nVoxels,nCoils,2);

% Estimate the best PD for each position a linear sulotion
% This makes it a nested bilinear problem
for kk=1:2

    %estime PD oh half coils
PD = zeros(nVoxels,1);
for ii=1:nVoxels
    use= find(useX(ii,:)==kk);
    PD(ii) = G(loc(ii),use)' \ M0(loc(ii),use)';
end


%calculate error on the other half
M0P = G(loc,:).*repmat( PD,1,nCoils);
err_ = M0_XV(useX~=kk) - M0P((useX~=kk));
err_X=[err_X err_];

end

err_X=err_X(:);



