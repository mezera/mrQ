function est=solveRatio(r,pMatrix)


Npar=size(pMatrix,2);
Ncoils=size(r,2);
Nvoxels=size(r,1);

% The big matrix has the polynomials from the  coils side by side.
% The first coil is always just the polynomail matrix
% The second coil is multiplied on every row by the ratio of the M0 values.
% The gain solutions satisfy 0 = rhs g
st_vox=1;
ed_vox=Nvoxels;

for i=1:Ncoils-1
    for j=i+1:Ncoils
% the i coil
             stP=1+(i-1)*Npar;
            edP=stP+Npar-1;
            rhs(st_vox:ed_vox,stP:edP)=[pMatrix];
            
            %the j coil
            stP=1+(j-1)*Npar;
            edP=stP+Npar-1;
            rhs(st_vox:ed_vox,stP:edP)= [diag(-r(:,i)./r(:,j))*pMatrix];
   
        % move to the next mesuraments
        st_vox=ed_vox+1;
        ed_vox=ed_vox+Nvoxels;
    end
 end
%solve
[U, ~] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);

