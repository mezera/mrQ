%%% Figure 4
%
% Illustrate the poyinomyal order of coils with phantom data
%
% AM/BW Vistaosft Team, 2013



%%  Make sure mrQ is on the path
addpath(genpath(fullfile(mrqRootPath)));
%% Generate example parameters for the coils from the phantom data

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 3;% Which box 
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.


%% get the polynomyails error


for nSamples=2:10 % The box is -nSamples:nSamples
for pOrder=1:3 %  polynomial order

phantomP(nSamples,pOrder) = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
end
end

%% get the polynomyails error

for nSamples=2:10 % The box is -nSamples:nSamples
for pOrder=1:3 %  polynomial order
    
    Voulume(nSamples,pOrder)=((nSamples*2+1)*2)^3; % the box we used are voulume with -nSamples:nSamples = (nSamples*2+1) voxel  a side.
    % the resultion of the phantom scan is 2mm. so multipal by 2. and ^3
    % for voulume
PE(nSamples,pOrder)=phantomP(nSamples,pOrder).percentError;
end
end


mrvNewGraphWin;
hold on


plot(Voulume(2:10,1),PE(2:10,1) ,'-k*');
plot(Voulume(2:10,2),PE(2:10,2) ,'-ko');
plot(Voulume(2:10,3),PE(2:10,3) ,'-ks');
legend('1st order' , '2nd order', '3rd order','Location','NorthWest')
xlabel(' Voulume mm^3');ylabel('pracent error');

%% 








