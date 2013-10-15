%% Fig_brainRepeatability
%
% Compares the same brain measured using 8 and 32 channel coil data
%
% To test the reliabilty of our method we scan the same subject in two
% different scans with 8 and 32 chanel coils. we compare the PD estimate
% in the two scans
%
% AM/BW Vistaosft Team, 2013

%% load images we need to put this data somewhere so it can be available
%(the library?)
%
% Directory on Windows.  There will be another for Mac.
% biacName = '\\red.stanford.edu\home\brian\BIAC\biac\biac2\wandell2\data'
%
%
% those are the Water fraction maps  (WF =~PD)  and the brain mask (BM)
% voulums o f the same subject scan in 32 and 8 chanels coils. the T1
% wighted images of each scan were register to wach other by ANTS sfotwere
% and the transform was applied on the BM and the WF volumes.

%% BW's computer (Windows)
thisFile8       = '/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/maps/WF_map.nii.gz';
thisFile32      = '/WMDevo/code/CompareMaps/AM_8_32/WarpWF_map.nii.gz';
thisFileBM32    = '/WMDevo/code/CompareMaps/AM_8_32/WarpbrainMask.nii.gz';
thisFileBM8     = '/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz';
thisFileCompare = '/WMDevo/code/CompareMaps/AM_8_32/WarpbrainMask.nii.gz';

fname = fullfile(biacName,thisFile8);
WF8ch = niftiRead(fname);

fname = fullfile(biacName,thisFile32);
WF32chWarp = niftiRead(fname);

fname = fullfile(biacName,thisFileBM8);
BM8ch =niftiRead(fname);

fname = fullfile(biacName,thisFileBM32);
BM32chWarp = niftiRead(fname);

%% AVIV
WF8ch      = niftiRead('/biac2/wandell2/data/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/maps/WF_map.nii.gz');
BM8ch      = readFileNifti('/biac2/wandell2/data/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
WF32chWarp = readFileNifti('/biac2/wandell2/data/WMDevo/code/CompareMaps/AM_8_32/WarpWF_map.nii.gz');
BM32chWarp = readFileNifti('/biac2/wandell2/data/WMDevo/code/CompareMaps/AM_8_32/WarpbrainMask.nii.gz');

%% Set up brain mask and calculate error score (R2)

% Take the brain mask as the intersection of the two masks
BM = logical(BM8ch.data) && logical(BM32chWarp.data);

% This region is cut off at the edge of scan box
BM(:,:,1:52)   = 0;
BM(:,:,144:end)= 0;

% The coefficient of determination (R^2)
CV = (calccod(WF8ch.data(BM),WF32chWarp.data(BM))/100).^2;

%% Plot a  2D histogram of the error

% We plot the MTV = 1 - WF the macromolecular tisuue volume fraction.
nBins = 155;
[n,x,y] = mrQ_hist2d(1-WF8ch.data(BM),1-WF32chWarp.data(BM),nBins);
maxN = ceil(max(n(:))/10)*10;

mrvNewGraphWin;
image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
identityLine(gca);
axis square xy; axis image;
xlim([0. 0.4]); ylim([0. 0.4])
ylabel('MTV  32ch' ,'FontSize',16);
xlabel('MTV  8ch' ,'FontSize',16);
set(gca,'FontSize',16)
title( [' R^2= ' num2str(CV)])
grid on

%% END