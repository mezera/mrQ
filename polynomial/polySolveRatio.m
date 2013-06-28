function  [eGains,  polyRatioMat ] = polySolveRatio(M0,pBasis)
% Solve for the gain parameters of the polynomial coefficient
%
%  eGains = polySolveRatio(r,pBasis)
%
% M0: A matrix whose columns are the M0 data at the (x,y,z)
%     positions (or x,y, or x). The ratios of these values in these columns
%     scale the pMatrix.  If there are P position and N coils this matrix
%     is P x N.
% pBasis:  The matrix of polynomial basis functions over positions
%
% eGains:        The estimated gains
% polyRatioMat:  The large matrix such that 0 = polyRatioMat * gains
%
% BW, AM Copyright Vistasoft Team, 2013

% Number of gain parameters
Npar    = size(pBasis,2);

% The data are stored in the columns of the matrix, r.
Ncoils  = size(M0,2);
Nvoxels = size(M0,1);

% The big matrix, polyRatioMat, has the polynomials from the coils in
% blocks.  These are typically separated by blocks of zeros. The first coil
% is always just the polynomial matrix, without any scaling The second coil
% is multiplied on every row by the ratio of the M0 values. The gain
% solutions satisfy 0 = M g
st_vox = 1;
ed_vox = Nvoxels;

for ii=1:Ncoils-1  % For each coil up to one before the list
    for jj=(ii+1):Ncoils  % Compare it to the next coil
        % the ith coil start and end positions
        stP = 1 + (ii-1)*Npar;
        edP = stP + Npar - 1;
        polyRatioMat(st_vox:ed_vox,stP:edP)= pBasis;
        
        %The new start and end positions
        stP = 1 + (jj-1)*Npar;
        edP = stP + Npar - 1;
        % Multiply the rows of pMatrix by the ratios of the M0 data from
        % the relevant coils
        polyRatioMat(st_vox:ed_vox,stP:edP)= diag(-M0(:,ii)./M0(:,jj))*pBasis;
        
        % move to the next group of measurements
        st_vox = ed_vox + 1;
        ed_vox = ed_vox + Nvoxels;
    end
end

%% The rows of polyRatioMat are ratios of M0 measurements
%  These measurements are from two coils at a voxel.
%  We could scale the row to reflect the SNR of that ratio.
%  For example, if one of the M0 values is small, that would reduce our
%  confidence in the ratio.
%
% Note, we have a lot of coil and voxels.  So we can get rid of quite a few
% (even zeroing them out) and still solve this equation.

%% I am concerned that the gain coefficients are so different in size
%
% The DC coefficient is 1000 times bigger than the other ones.  
% In fact, they decline precipitously across order.  So, maybe we should be
% doing something to balance these better.
%% We could re-write the M0 matrix to represent virtual coils
%
%  Each column of M0 is a real coil.  If we analyze the 

% Solve for the eigenvectors of polyRatioMat.  The first one has the
% smallest singular value and is thus the best estimate.  We scale the gain
% vector so that the first entry (the constant) is 1.
[U, ~] = eig(polyRatioMat'*polyRatioMat);
eGains = U(:,1)/U(1,1);

% R=r(:,ii)./r(:,jj); F=find(R<2);polyRatioMat1=polyRatioMat;
% polyRatioMat1(F,:)=0;
% [U, ~] = eig(polyRatioMat1'*polyRatioMat1);
% est1 = U(:,1)/U(1,1);

% %is a different way to solve the matrix make a different.?
% lhs=-polyRatioMat(:,1);
% rhs=polyRatioMat(:,2:end);
% est1 = rhs\lhs;
% no it's not up to a constant 

return
