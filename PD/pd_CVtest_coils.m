function Fit=pd_CVtest_coils(PD,pBasis,M0_v,plotFlag,TruePar)

nCoils=size(M0_v,2);

for ii=1:nCoils
    G(mask,ii)  = M0_v(mask,ii) ./ PD(:);         % Raw estimate
    g(:,ii) = pBasis(:,:) \ G(:,ii);  % Polynomial approximation
    
end





M0 = G.*repmat( PD(:),1,nCoils);
%calculate the fit error
Fit.M0prederr=M0-M0_v;
Fit.Meanerr=mean(abs(M0prederr));

fprintf(' the mean abs error fit  % 0.f  \n :',Fit.Meanerr);

Fit.PD=G;
Fit.M0=M0;

%plot
if plotFlag==1
    
    mrvNewGraphWin;plot(Fit.M0prederr,'.');
    
    
    
    mrvNewGraphWin;
    
    if notDefined('TruePar')
        
    else
        Par=TruePar;
        str='True gains';
    end
    CoefNorm=Par(1)./g(1);
    
    % Plot the starting point
    figure(figH);
    set(figH,'Name', ['Loop ' num2str(k)]);
    subplot(1,2,1); plot(g(:).*CoefNorm,Par(:),'.')
    identityLine(gca);
    xlabel('Estimated gains'); ylabel(str);
    
    subplot(1,2,2); plot(PD,'.');set(gca,'ylim',[min(PD(:))*0.9 max(PD(:))*1.1]);
    ylabel('Estimated PD');xlabel('voxel');
    
end
if ~notDefined('TruePD')
    Fit.RMSE = sqrt(mean(  (g(:).*CoefNorm-Par(:)  ).^2));
    fprintf(' the coil coefficents RMSE % 0.f  \n :',  Fit.RMSE);
    
end

