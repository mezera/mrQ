
coilList = [1:3];
        Par=OutPut.params(:,coilList); Par=Par./Par(1);
Lamda=0.1;
tryagain=1;

% the right stat good guess
%PD=ones(size(OutPut.M0S_v(:,1)),1);

%sum of sqr of M0 as a starting point 
PD=   sqrt(sum(OutPut.M0_v(:,coilList).^2,2)  )   ;


 clear g0
 for ii=1:length(coilList)
     G=OutPut.M0_v(:,coilList(i))./PDn;
g0(:,ii)= OutPut.pBasis \ G;
 end
%PD=randn(size(OutPut.M0S_v(:,1)),1);

PD=PD./mean(PD);
k=0;

mrvNewGraphWin(num2str(k)); subplot(1,2,1);(plot(g0(:),Par(:),'.'))
        identityLine(gca);
        xlabel('Estimated gains'); ylabel('True gains');
        subplot(1,2,2);plot(PD,'.')

while tryagain==1
    g=RidgeRegressCoilfit(PD,Lamda,OutPut.M0_v(:,coilList),OutPut.pBasis);
    
    G=OutPut.pBasis *g;
    PDn=OutPut.M0_v(:,coilList) .\G;
    PDn=mean(PDn,2);
    PDn=PDn./mean(PDn(:));
    k=k+1;
    if std(PD-PDn)<0.01
        tryagain=0;
    else
        mrvNewGraphWin(num2str(k)); subplot(1,2,1);(plot(g(:),Par(:),'.'))
        identityLine(gca);
        xlabel('Estimated gains'); ylabel('True gains');
        subplot(1,2,2);plot(PDn,'.')
        
        PD=PDn;
        
    end
    if k==100
        tryagain=0;
    end
end


PDFinal=(PD+PDn)/2;
PDFinal=reshape(PDFinal,OutPut.SZ(1:3));
showMontage(PDFinal)
