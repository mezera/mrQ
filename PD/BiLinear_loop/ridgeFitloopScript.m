
coilList = [1:8];
Par=OutPut.params(:,coilList); Par=Par./Par(1);
Lamda=0.1;
tryagain=1;
M0_v=OutPut.M0_v(:,coilList);

% simlatedPD=OutPut.M0S_v(:,4);
% simlatedPD=simlatedPD./mean(simlatedPD);
% for ii=1:length(coilList)
%     M0_v(:,ii)=M0_v(:,ii).*simlatedPD;
% end
%     PDtrue=reshape(simlatedPD,OutPut.SZ(1:3));
%     showMontage(PDtrue);

% the right stat good guess
%PD=ones(size(OutPut.M0S_v(:,1)),1);

%sum of sqr of M0 as a starting point
PD=   sqrt(sum(M0_v.^2,2)  )   ;
PD=PD./mean(PD);


clear g0 G
for ii=1:length(coilList)
    G(:,ii)=M0_v(:,ii)./PD;
    g0(:,ii)= OutPut.pBasis \ G(:,ii);
end
%PD=randn(size(OutPut.M0S_v(:,1)),1);

k=0;

mrvNewGraphWin(num2str(k));
subplot(1,2,1);plot(g0(:),Par(:),'.')
identityLine(gca);
xlabel('Estimated gains'); ylabel('True gains');
subplot(1,2,2);plot(PD,'.')
ylabel('Estimated PD');xlabel('voxel');

while tryagain==1
    g=RidgeRegressCoilfit(PD,Lamda,M0_v,OutPut.pBasis);
    
    Gn=OutPut.pBasis *g;
    PDn=M0_v .\Gn;
    PDn=mean(PDn,2);
    PDn=PDn./mean(PDn(:));
    k=k+1;
    %        if std(PD-PDn)<0.01
    if std(G-Gn)<0.004
        tryagain=0;
    else
        PD=PDn;
        G=Gn;
        mrvNewGraphWin(num2str(k));
        subplot(1,2,1);plot(g(:),Par(:),'.')
        identityLine(gca);
        xlabel('Estimated gains'); ylabel('True gains');
        subplot(1,2,2);plot(PDn,'.')
        ylabel('Estimated PD');xlabel('voxel');
        
        
    end
    if k==100
        tryagain=0;
    end
end


PDFinal=(PD+PDn)/2;
PDFinal=reshape(PDFinal,OutPut.SZ(1:3));
showMontage(PDFinal)
