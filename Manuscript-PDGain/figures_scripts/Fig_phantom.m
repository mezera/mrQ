%%% Figure 4
%
% Illustrate the polynomial order of coils with phantom data.  In this case
% we assume 
%
% AM/BW Vistaosft Team, 2013

%%  Make sure mrQ is on the path
addpath(genpath(fullfile(mrqRootPath)));

%% Set parameters for the coils from the phantom data

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 3;% Which box
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.


%% Get the polynomial error

% We load the phantom data.
% We fit the polynomial error for a range of box sizes and polynomial
% orders.  
% In the phantom test case we know that PD is constant.
for nSamples=2:10         % The box is -nSamples:nSamples
    for pOrder=1:3        %  polynomial order
        phantomP(nSamples,pOrder) = ...
            pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
            noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
    end
end

%% Extract key values from the structure containing the polynomial fits

% Number of samples is the box size
% pOrder is the polynomial order
PE = zeros(10,30);
Volume = zeros(10,3);
for nSamples=2:10 % The box is -nSamples:nSamples
    for pOrder=1:3 %  polynomial order
        
        % The box we use is a volume with -nSamples:nSamples on a side.
        % The resolution of the phantom scan voxel is 2mm. 
        % So multply by 2 and ^3 for volume.
        % The volume in mm3 is 
        Volume(nSamples,pOrder) = ((nSamples*2+1)*2)^3; 
        PE(nSamples,pOrder)     =  phantomP(nSamples,pOrder).percentError;
    end
end

% Show the fit
mrvNewGraphWin;
hold on
plot(Volume(2:10,1),PE(2:10,1) ,'-k*', 'MarkerSize',10);
plot(Volume(2:10,2),PE(2:10,2) ,'-ko', 'MarkerSize',10);
plot(Volume(2:10,3),PE(2:10,3) ,'-ks', 'MarkerSize',10);
legend('1st order' , '2nd order', '3rd order','Location','NorthWest')
xlabel(' Volume mm^3','FontSize',16);ylabel('percent error','FontSize',16);
set(gca,'xscale','log','yscale','linear','FontSize',16)
axis image; axis square

% So, we are using 3rd order polynmials with nSamples = 7, which is 1.4 cm
% on a side.  Goldilocks.
%
% Question - should we run some cross-validation? Maybe, maybe not.
%
% The idea: We hold out, say half the data and fit a polynomial to one half
% and show the error on the other half.  We choose the volume that has the
% smallest cross-validation error fits.  We then fit that and move on.
%

%% End

