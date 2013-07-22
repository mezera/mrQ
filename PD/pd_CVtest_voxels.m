function Fit=pd_CVtest_voxels(g,pBasis,M0_v,plotFlag,TruePar,TruePD,FigNum)
% Make X-validation values and plots
%
%
%
%
% AM VISTASOFT Team, 2013
if notDefined('plotFlag'), plotFlag = 0; end

nCoils=size(M0_v,2);

% Estimate PD from the coil gain and M0 data
[PD, G] = pdEstimate(M0_v, pBasis, g);

% Normalize to first PD entry of 1
G  = G  .* PD(1);
PD = PD ./ PD(1);

%Get the M0 prediction
M0 = G.*repmat( PD(:),1,nCoils);

%Calculate the fit error voxel by voxel
Fit.M0prederr = M0 - M0_v;

% This is average cross each coil
Fit.Meanerr = mean(abs(Fit.M0prederr));
fprintf(' the M0 mean abs error fit  % 0.3f  \n :',Fit.Meanerr);

if ~notDefined('TruePD')
    Fit.RMSE = sqrt(mean(  (TruePD(:)./mean(TruePD(:))-PD(:)./mean(PD(:))   ).^2));
    fprintf(' the PD  RMSE % 0.3f  \n :',  Fit.RMSE);
end

Fit.PD = PD;
Fit.M0 = M0;

if ~notDefined('FigNum')
    fprintf('Fig number % 0.f  \n :',FigNum);
end

% Plot
if plotFlag==1
    if notDefined('FigNum'), FigNum=1; end
    figH = mrvNewGraphWin([],'tall');
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
        subplot(3,1,1); plot(g(:).*CoefNorm,Par(:),'.')
        identityLine(gca);
        xlabel('Estimated gains'); ylabel(str);
        
        subplot(3,1,2); plot(PD(:),TruePD(:),'.');
        set(gca,'ylim',[min(PD(:))*0.9 max(PD(:))*1.1]);
        xlabel('Estimated PD');ylabel('true PD');
        title(['the PD RMSE: ' num2str(Fit.RMSE)])
        
        subplot(3,1,3);
        plot(Fit.M0prederr./M0_v,'.');
        title(['the M0 mean abs error : ' num2str(Fit.Meanerr)])
        
    end
end