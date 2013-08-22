%% Compute the error for different polynomial orders and box volumes
%
% We use this script to write out the error and choose the order and size
% for other applications.
%
% The data used to assess are from the phantom.
% These are 32 coils at the CNI 3T.
%
% Conclusion:  The 3rd order with 10 samples (2 cm) is fine.
%              The 2nd order is good up to about 5 samples (1 cm)
%              The 1st order is rarely OK, and often poor
%
% AM/BW Copyright VISTASOFT, 2013

% The volume is =nSamples:nSamples.
% Each voxel in the data is 2 mm isotropic


%% Loop on sizes and polynomial orders

nCoils     = 32;
nDims      = 3;
noiseFloor = 500;
sampLocation = 3;
maxSample    = 12;
maxPolyOrder = 3;
percentError = zeros(maxSample,maxPolyOrder); 

% We need at least 5 samples to go up to 3rd order.
% We could check 1st and 2nd order with fewer samples.
for nSamples=2:maxSample
    % Polynomial order
    for pOrder = 1:maxPolyOrder
        [~,~, ~,~,percentError(nSamples,pOrder)] = ...
            pdPolyPhantomOrder(nSamples,nCoils,nDims,pOrder,noiseFloor,sampLocation);
    end
end



%% Plot the results
mrvNewGraphWin;
imagesc(1:maxPolyOrder, 1 + (1:maxSample)*2,log10(percentError(2:end,:)*100))
colormap(hsv)
xlabel('Polynomial order')
ylabel('Number of samples per dimension')
zlabel('Percent error')
cb = colorbar;
title('Log10 percent error')


% yTicks = (1:5:2*maxSample);
% set(gca,'yTick',(1:5:2*maxSample),'yTickLabel',fliplr(yTicks))


%%
mrvNewGraphWin;
semilogy(1 + (2:maxSample)*2, percentError(2:end,:)*100,'-*')
legend('pOrder_1','pOrder_2','pOrder_3')
xlabel('nSamples')
ylabel('Log 10 percent error')
grid on
ln = line([0 30],[1 1]); set(ln,'Color', [ 0 0 0])


%% Save the file
name = (fullfile(mrqRootPath,'PD','figures', 'percentError_nova32Ch'));
save(name,'VarEx');
