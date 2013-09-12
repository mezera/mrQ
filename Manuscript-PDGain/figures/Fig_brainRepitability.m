% figure 6
%
% To test the relaibilty of our method we scan the same subject in two
% different scan with 8 and 32 chanel coils. we compare the PD estimation
% in the two scans
%
%
% AM/BW Vistaosft Team, 2013
%%
% load images we need to put this data somewhere so it can be avilable 
%(the librerey?)

% those are the Water fraction maps  (WF =~PD)  and the brain mask (BM) voulums o f the same subject scan in 32 and 8 chanels coils. the T1 wighted images of each scan were
% register to wach other by ANTS sfotwere and the transfrom was applay on the
% BM and the WF volumes.

WF8ch=readFileNifti('/biac2/wandell2/data/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/maps/WF_map.nii.gz')
WF32chWarp=readFileNifti('/biac2/wandell2/data/WMDevo/code/CompareMaps/AM_8_32/WarpWF_map.nii.gz');

BM8ch=readFileNifti('/biac2/wandell2/data/WMDevo/adult/108_AM/QuantitativeImaging/20111020_1297_8ch_1mm3/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
BM32chWarp=readFileNifti('/biac2/wandell2/data/WMDevo/code/CompareMaps/AM_8_32/WarpbrainMask.nii.gz');

BM=logical(BM8ch.data) & logical(BM32chWarp.data);

% this region are cuted or adge of scan box
BM(:,:,1:52)=0;
BM(:,:,144:end)=0;

%  the coefficient of determination (R^2) 
CV=(calccod(WF8ch.data(BM),WF32chWarp.data(BM))/100).^2;

%% plot a  2D histogram 
% we plot the 1-WF= MTV the macromulecular tisuue voulume fraction.

nBins = 155;
[n,x,y] = mrQ_hist2d(1-WF8ch.data(BM),1-WF32chWarp.data(BM),nBins);
maxN = ceil(max(n(:))/10)*10;


mrvNewGraphWin;

image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
 identityLine(gca); 
axis square xy;
axis image; 
xlim([0. 0.4]); ylim([0. 0.4])
           ylabel('  MTV  32ch' ,'FontSize',16);

           xlabel('   MTV  8ch' ,'FontSize',16);
set(gca,'FontSize',16)

title( [' R^2= ' num2str(CV)])

% end