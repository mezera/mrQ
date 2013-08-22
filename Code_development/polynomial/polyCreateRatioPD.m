function  [polyRatioMat, M0pairs] = polyCreateRatioPD(M0,pBasis)
% Solve for the gain parameters of the polynomial coefficient
%
%  [polyRatioMat, M0pairs] = polyCreateRatioPD(r,pBasis)
%
% M0: A matrix whose columns are the M0 data at the (x,y,z)
%     positions (or x,y, or x). The ratios of these values in these columns
%     scale the pMatrix.  If there are P position and N coils this matrix
%     is P x N.
% pBasis:  The matrix of polynomial basis functions over positions
%
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
ed_vox = Nvoxels*2;
M0pairs = [];

for ii=1:Ncoils-1  % For each coil up to one before the list
    for jj=(ii+1):Ncoils  % Compare it to the next coil
        
        % the iith coil start and end positions
        stP = 1 + (ii-1)*Npar;
        edP = stP + Npar - 1;
        polyRatioMat(st_vox:ed_vox-Nvoxels,stP:edP)= pBasis;
        
        polyRatioMat(ed_vox-Nvoxels+1:ed_vox ,stP:edP)=diag(M0(:,jj))*pBasis;

         % the jjth coil start and end positions
        stP = 1 + (jj-1)*Npar;
        edP = stP + Npar - 1;
        
        % Multiply the rows of pMatrix by the ratios of the M0 data from
        % the relevant coils
        polyRatioMat(st_vox:ed_vox-Nvoxels,stP:edP)= diag(-M0(:,ii)./M0(:,jj))*pBasis;
        M0pairs = [M0pairs; [M0(:,ii), M0(:,jj)]]; %#ok<AGROW>
        polyRatioMat(ed_vox-Nvoxels+1:ed_vox ,stP:edP)=-diag(M0(:,ii))*(pBasis);

        % move to the next group of measurements
        st_vox = ed_vox + 1;
        ed_vox = ed_vox + Nvoxels*2;
    end
end


return
