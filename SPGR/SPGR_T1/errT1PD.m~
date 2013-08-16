function err =errT1PD(x,flipAngles,tr,S,Gain,B1,lsq,SD)
%
% err =errT1PD(x,flipAngles,tr,S,Gain,B1,lsq,SD)
%
% estimate the fit of x(1) PD and x(2) T1 to fit the SPGR T1 data with
% different flip angles
%
%InPut
% x the fitted parameters
% flipAngles the scansflipAngles 
% tr the scans tr
% S - the mesured SPGRs images
% Gain the coil gain (can be also set to be one and fitted later
% B1 the exsite inhomogenity (the error in nominal flipangles
%lsq the kind of error calculation
% SD a way to normalize the error by the data std.
%
%OutPut
%err -the error between the estimation and the data.
     
  
     if (~exist('lsq','var')||isempty(lsq)) 
    lsq=1;
end;
if lsq~=0, lsq=1;end
    if (~exist('SD','var')||isempty(SD)),
    SD=1;
end;
M0=x(1).*Gain;
fa=flipAngles.*B1;
fa = fa./180.*pi;
% the SPGR eqation 
Sc =M0.*(1-exp(-tr./x(2))).*sin(fa)./(1-exp(-tr./x(2)).*cos(fa));


if lsq==1
    err=1/SD.*(S-Sc);
    err=sqrt(abs(err)); %let fit the median and not the mean that will give less whiat for outlayers

else,

err=1/SD.*sum(sum(abs ((S-Sc)./S) ));

 end;
