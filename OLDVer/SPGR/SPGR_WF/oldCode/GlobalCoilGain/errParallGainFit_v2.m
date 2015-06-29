function [err]=errParallGainFit_v2(x,box,Poly,coefdat,coefdatT1,mask2)
%
%[err]=errlocalGainUC_v2(x,box,Poly,coefdat,use,coils)
%
% # estimation of the coils bias. we try to fit the different bias so:
% 1. all coil will have the same pd -->  PD=Coil_data/Coil_Bias
% 2. the different coils should not gain corrlation after removing the pd
% from all of them. becouse grater coralation mean adding bias part that
% will exsist in all coils and therefor pd estimations
%
% InPuts:
% x     - are the polynomials coeficient for each coil 
% his is a 2D matrix ( coeficient X coils )  
%
% box   -  is the coils data we fit to (2d (coils X data) (box is about 20^3x10).
%
% Poly  - are the polynomials we multipal the coeficient with
%
% coefdat -  coralation coeficient of the input coil by coil data ( in box).
%
% use - the location of the usiable data along the box (this is brain mask and clean for  out layers).
%
%OutPut
%err -the error in the estimation for unique PD with no added bias

%%
%arbitrary Penalty parameter
correlationPenalty = 5;

% The gain for each of the coils, estimated by the parameter x
Gain = Poly*x(:,:);
%Gain=Gain(mask,:);
if find(size(Gain)~=size(box))
    size(Gain)
    size(box)
    err=ones(1,length(Poly))+1e10;
else    
% TV=M0/Gain
box=box./Gain;

 %err0=mean(Gain(:,1))-MG1;
 
 err0=   nanmean(box(:))-1;

%%  calculate the PD
%Gain =reshape(Gain(use),[],coils);
%this is the predicted brain in each coil (row)
%  the reshape is faster way to do that --> the gain in the useable voxels
%    for i=1:10; tt(:,i)=Gain(use(:,i),i);end
%box=box./Gain(mask,:)
%Val = box./Gain;
%Val
%Val(~M0mask)=0;
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

%%
%%%the error in fit
%err1=(abs(box - repmat(nanmean(box,2),1,coils)) )./repmat(nanmean(box,2),1,coils);
err=(nanstd(box,[],2));

%err1=( abs(Val - repmat(mean(Val,2),1,coils)) )./repmat(mean(Val,2),1,coils);
%err2=std(err1,[],1);


%err1=sum(err1,2);

%% we don't use this for now
%we prefer a constant error and not local fits
%err1=err1*(1+sum(err2));

%%
% %, turn your attention to the correlation (overlap) between the coil
% gain functions.  We compute the corrcoefs and and take out the lower
% triangular (non-redundant) part of these.  -1 means take everything below
% the diagonal

%Gain(~M0mask)=0;
Gain(isnan(box))=0;
coefG =  tril(corrcoef(Gain),-1);

Gain(~mask2,:)=0;
coefGT1 =  tril(corrcoef(Gain),-1);



% Outside of this routine we already calculated
%   coefdat =tril(corrcoef(box),-1);
% Here we compare the two correlation coefficients
% The correlation coefficients of the data (coefdat) should be larger than
% the corr coef of the gain (coefG).
err2 = (coefdat(coefdat~=0)-coefG(coefdat~=0))./abs(coefdat(coefdat~=0)); 
err22 = (coefdatT1(coefdatT1~=0)-coefGT1(coefdatT1~=0))./abs(coefdatT1(coefdatT1~=0)); 

% If the difference is positive, there is no penalty.
err2(err2>0)=0;
err22(err22>0)=0;

% The remaining terms have a penalty, so we add them up and multiply by 5.
% That value is arbitrary.
%err2 = sum(abs(err2))*correlationPenalty;
%err22 = sum(abs(err22))*correlationPenalty;

err2 = (abs(err2))*correlationPenalty;
err22 = (abs(err22))*correlationPenalty;




% Combine the two errors by multiplying them.
% Could be err1 + correlationPenalty*err2  
%err = err1*(1+err2);
%err = err1*(1+err22);

%% check for crazy fit
%let not alow minos brain!!!
err3=0;
if find(box<0) ; err3=mean(box(box<0));


err=err*(1+abs(err3));end
%% 
%gain shouldn't go under 1 there is no such thing

%we checak if the coil gain fit goes wild. Gain can't be less then one we
%don't lose signal ....
err4=0;
if find(Gain<1 & Gain~=0) ; err4=1-mean(Gain(Gain<1 & Gain~=0));
err=err*(1+abs(err4));
end

%we try to keep the mean of the gain constat 
%err=err*(1+abs(err0));
err=[err' err2'*20 err22'*20 err0*100];
%% we don't use this for now

%let not allow crazy local fit
% SDF=std(Val,2);
% if SDF>SDD*10
%     err=err*(1+abs(SDF/SDD));
% end;

%%
% Replace the infinite values with a large value to help later computations
err(isinf(err))=1e10;
err=double(err);
end
return

