function Fit=pd_CVtest_coils(PD,pBasis,M0_v,plotFlag,TruePar,FigNum)

nCoils=size(M0_v,2);

for ii=1:nCoils
    G(:,ii)  = M0_v(:,ii) ./ PD(:);         % Raw estimate
    g(:,ii) = pBasis(:,:) \ G(:,ii);  % Polynomial approximation
    
end


M0 = G.*repmat( PD(:),1,nCoils);
%calculate the fit error
Fit.M0prederr=M0-M0_v;
Fit.Meanerr=mean(abs(Fit.M0prederr));

fprintf(' the mean abs error fit  % 0.f  \n :',Fit.Meanerr);

 if ~notDefined('TruePar')
     CoefNorm=TruePar(1)./g(1);
        Fit.RMSE = sqrt(mean(  (  g(:).*CoefNorm-TruePar(:)   ).^2) );
    fprintf(' the coil coefficent  RMSE % 0.3f  \n :',  Fit.RMSE);
    
end


Fit.PD=G;
Fit.M0=M0;

%plot
if plotFlag==1
    if notDefined('FigNum')
        FigNum=1;
    end
    figH = mrvNewGraphWin;
    set(figH,'Name',   num2str(FigNum));
    
    if notDefined('TruePar')
        
        plot(Fit.M0prederr./M0_v,'.');
        xlabel('voxels');ylabel('M0 pracent error')
        title(['the M0 mean abs error : ' num2str(Fit.Meanerr)])
    else
        
        Par=TruePar;
        str='True gains';
        
        CoefNorm=Par(1)./g(1);
        
        % Plot the starting point
        figure(figH);
        subplot(1,2,1); plot(g(:).*CoefNorm,Par(:),'.')
        identityLine(gca);
        xlabel('Estimated gains'); ylabel(str);
       title([' the coil coeficents error ' num2str(Fit.RMSE) ])
        
        subplot(1,2,2);
        plot(Fit.M0prederr./M0_v,'.');
        title(['the M0 mean abs error : ' num2str(Fit.Meanerr)])
        
        
    end
end



