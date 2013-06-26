function[ Res]   = fitRatioandPlotPD(coilList,M0_v,M0,pMatrix,params, name,Figs)

if notDefined('name')
    name='';
end
if notDefined('Figs')
    Figs=0;
end

[Res.est Res.polyRatioMat] = polySolveRatio(M0_v(:,coilList),pMatrix);

 [Res.CoilCoefErr Res.PDerr  Res.estMatrix  Res.ParMatrix Res.G Res.M00 Res.PD ...
    Res.PDspaceErr]=polyRatioErr(Res.est,params(:,coilList),pMatrix);

Res.R=M0(:,:,:,coilList(1))./M0(:,:,:,coilList(2));
MinPD=min(Res.PD(:));
MaxPD=max(Res.PD(:));

if Figs==1
showMontage(Res.PD,[],[],[],[],1);colormap hot; title(['PD' name]);caxis([MinPD MaxPD]);
 showMontage(Res.R,[],[],[],[],2);colormap hot;title(['Ratio' name])
 showMontage(M0(:,:,:,coilList(1)),[],[],[],[],3);colormap hot;title(['C1' name])
showMontage(M0(:,:,:,coilList(2)),[],[],[],[],4);colormap hot;title(['C2' name])
 showMontage(Res.PDspaceErr,[],[],[],[],5);colormap hot;title(['PD error' name])
figure(6);plot(Res.R(:),Res.PD(:),'*'); xlabel('R');ylabel(['PD' name]);
end



