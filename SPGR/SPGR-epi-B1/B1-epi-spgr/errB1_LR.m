function err =errB1_LR(x,flipAngles,tr,S,T1,f,mask,ratios)
%   
%err =errB1_LR(x,flipAngles,tr,S,T1,f,mask,ratios)
% estimate the fit B1 to voxelsthat have a known T1 ( from SEIR) and SPGR (S) measres.
%the SPGR are data mesure with different flip angles. x is a single B1
%value that is assign to several voxel (this will generate smothdness) for
%overlap solusions, the different voxel are wighted in space by f.
%
%Arrguments
% x the fitted parameter B1 paraeter: 
% flipAngles the scansflipAngles 
% tr the scans tr
% S - the mesured SPGRs images   (voxel,flipangle,tr)
% T1  the T1 value calculated from SEIR  (voxel,1)
% f is waithed function 

%outPut
%err -the error between the estimation and the data.
     
B1 = x;


% the SPGR eqation S  Sc=S./M0 
for ii=1:length(tr)
    fa=flipAngles(ii).*B1;
fa = fa./180.*pi;

Sc(:,ii) =(1-exp(-tr(ii)./T1(:))).*sin(fa(:))./(1-exp(-tr(ii)./T1(:)).*cos(fa));

end

measure=S(mask,ratios(:,1))./S(mask,ratios(:,2));
model=Sc(mask,ratios(:,1))./Sc(mask,ratios(:,2));

err=((measure-model));
err=err.*f(mask,:); % waiting by the filter (local invierment)

err=sqrt(abs(err)); %let fit the median and not the mean that will give less whiat for outlayers









