function [err]=errlocalGainUC_v1(x,box,Poly,coefdat,use,coils)
% x are the polynomials coeficient. size(x) =(10,number of coils). 10 paramters define is 2ed order 3D poly.  
% box is 2 d matrix the input mesuraments the 1st dimention is the mesurament in the (brain mask --use) voxels. the 2ed dimantion is the 
% different coils number. box is about 20^3x10.
% Poly are the polynomials
% coefdat coralation coefs of the data.
% use location of the usiable data along box (this is brain mask and no out layers).
%fitting the local Gain, we don't allow coralation between the coils% but we

% static?
correlationPenalty = 5;

% The gain for each of the coils, estimated by the parameter x
Gain = Poly*x(:,:)';
Gain =reshape(Gain(use),[],coils);
%this is the predicted brain in each coil (row)
%  the reshape is faster way to do that --> the gain in the useable voxels
%    for i=1:10; tt(:,i)=Gain(use(:,i),i);end
Val = box./Gain;
% We want the mean to stay invariant, only minimizing on the total
% variation, not the mean
%Val = Val./mean(Val(:)); %we normalize by the mean 

% This is the std of the estimate of each brain voxel.  We want the
% coefficients find a solution that minimizes the different brain estimates
% across the coils. So this std is across the coil dimension.  We get a std
% for every box entry (brain voxel).
%err1 = std(Val,[],2);
% faster with norm(.,1) - this means sum(abs(val - mean(val)))
%
%err1=sum((abs(Val - mean(Val(:)))),2);
err1=sum(((abs(Val - repmat(mean(Val,2),1,coils)) )./repmat(mean(Val,2),1,coils)),2);

% %, turn your attention to the correlation (overlap) between the coil
% gain functions.  We compute the corrcoefs and and take out the lower
% triangular (non-redundant) part of these.  -1 means take everything below
% the diagonal
coefG =  tril(corrcoef(Gain),-1);

% Outside of this routine we already calculated
%   coefdat =tril(corrcoef(box),-1);
% Here we compare the two correlation coefficients
% The correlation coefficients of the data (coefdat) should be larger than
% the corr coef of the gain (coefG).
err2 = (coefdat(coefdat~=0)-coefG(coefdat~=0))./abs(coefdat(coefdat~=0)); 

% If the difference is positive, there is no penalty.
err2(err2>0)=0;

% The remaining terms have a penalty, so we add them up and multiply by 5.
% That value is arbitrary.
err2 = sum(abs(err2))*correlationPenalty;

% Combine the two errors by multiplying them.
% Could be err1 + correlationPenalty*err2  
err = err1*(1+err2);
err3=0;
if find(Val<0) ; err3=mean(Val(Val<0));end


err=err*(1+abs(err3));
% Replace the infinite values with a large value to help later computations
err(isinf(err))=1e10;


%and let not alow minos brain!!!
return







err0=((abs(Val - repmat(mean(Val,2),1,coils)) )./repmat(mean(Val,2),1,coils));
% 
Gain1 = Poly*x(:,:)';
 box1=zeros([27 27 25 coils]);
 val1=zeros([27 27 25 coils]);
 G=zeros([27 27 25 coils]);
  val2=zeros([27 27 25 coils]);
for i=1:coils,
    
    tmp=zeros([27 27 25]);
    tmp(use(:,1))=Val(:,i);
val1(:,:,:,i)=tmp;
  tmp(use(:,1))=box(:,i);
box1(:,:,:,i)=tmp;
tmp=reshape(Gain1(:,i),size(tmp));
G(:,:,:,i)=tmp;
val2(:,:,:,i)=box1(:,:,:,i)./G(:,:,:,i);
end


%%
% a simultion for the coraltion reduction argument
% a1 = posrect(3+randn(500,1000));
%  a2 = posrect(3+randn(500,1000));
%  a3 = posrect(3+randn(500,1000));
%  a1a2 = calccorrelation(a1,a2,2);
%  a1a2new = calccorrelation(a1.*a3,a2.*a3,2);
%  figure;scatter(a1a2,a1a2new,'r.');
%  xlabel('original correlation (r)');
%  ylabel('correlation after corruption by a3 (r)');
