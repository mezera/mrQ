  

options = optimset('Display','iter','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun', 1e-10,'TolX', 1e-10);
  
  
  
  x0=BLSim.g;
[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par)  errFitNestBiLinear(par,M0SN,OutPut.pBasis,nVoxels,nCoilsS)...
         ,double(x0),[],[],options);
G = OutPut.pBasis*res1;

PD = zeros(nVoxels,1);
for ii=1:nVoxels
    PD(ii) = G(ii,:)' \ M0SN(ii,:)';
end
PDfit = reshape(PD,OutPut.SZ(1:3));
showMontage(PDfit);

showMontage(PDsim./mean(PDsim(:))-PDfit./mean(PDfit(:))  );
sum(abs(PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))))

RMSE = sqrt(mean(  (PDsim(:)./mean(PDsim(:))-PDfit(:)./mean(PDfit(:))   ).^2))
title(['the percent error    RMSE = '   num2str(RMSE)] )
